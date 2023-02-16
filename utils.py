#!/usr/bin/env python3

def extract_clusters(assoc_grp_df, window):
    
    import pandas as pd
    
    cands = []
    cand = []
    bp_list = assoc_grp_df.BP.tolist()
    bp_length = len(bp_list)
    for idx, pos in enumerate(bp_list):
        
        if idx + 1 == bp_length:
            if len(cand) > 1:
                cands.append(tuple(cand))
            break
            
        if bp_list[idx + 1] - bp_list[idx] > window:
            if len(cand) > 1:
                cands.append(tuple(cand))
                cand.clear()
            continue
        
        if idx not in cand:
            cand.append(idx)
        cand.append(idx+1)
    
    return cands

def efficient_query(assoc_df, chrom, start, end):
    
    import pandas as pd
    
    search_start = assoc_df["CHR"].searchsorted(int(chrom))
    search_end = assoc_df["CHR"].searchsorted(int(chrom) + 1)
    
    by_chrom = assoc_df.iloc[search_start: search_end, :]
    pos_start = by_chrom.BP.searchsorted(start) + by_chrom.index[0]
    pos_end = by_chrom.BP.searchsorted(end) + by_chrom.index[0]
    
    return assoc_df.iloc[pos_start:pos_end+1, :]

def second_pass_singletons(snp_df, outdir, window=5000, p_thresh=1e-3):
    '''
    This function implements second-pass filter on SNPs with p < 0.05.
    Signal aggregates with singletons (an isolated variant with very 
    high significance) are removed. Significant variants and clusters
    they lie in are reported in snps.id, which is a list of such var-
    iants. By default, singletons are found with a 5kb window. p_thr-
    esh is a user-defined parameter that controls the level of signi-
    ficance for a variant to be considered as significant.
    
    Returns:
    --------
    clusters - list of Cluster objects
    '''
    
    import os
    import pandas as pd
    import variants
    import logging
    
    logger = logging.getLogger("root") 
    
    clusters, singletons, sig_vars_cnt, sig_cluster_cnt = [], 0, 0, 0
    fout = open(os.path.join(outdir, "snps.id"), "w")
    for cluster in extract_clusters(snp_df, window):
        var_lst = []
        for idx, var in snp_df.iloc[list(cluster), :].iterrows():
            var_lst.append(variants.Variant(var.BP, var.P, var.CHR, var.SNP))
        cluster_cand = variants.Cluster(var_lst)
        if cluster_cand.get_outliers() == 1:
            singletons += 1
#             logger.debug("Identified one cluster classified as singleton. ")
        else:
            clusters.append(cluster_cand)
            if cluster_cand.minp < p_thresh:
                sig_cluster_cnt += 1
                for var in cluster_cand.get_significant_vars(p_thresh):
                    sig_vars_cnt += 1
                    fout.write(var.snpid)
                    fout.write("\n")
    fout.close()
    logger.debug(f"Identified {singletons} singletons misclassified as clusters. ")
    logger.info(f"Number of analytic clusters passing first singleton filter: {len(clusters)} ")
    logger.info(f"Number of significant variants (p < {p_thresh}): {sig_vars_cnt}")
    logger.info(f"Number of clusteres with sig variants (p < {p_thresh}): {sig_cluster_cnt}")
    logger.info(f"SNP ID list written to {os.path.join(outdir, 'snps.id')} ")
    
    return clusters

def third_pass_ld(ld_fp, all_snp_df, clusters, p_thresh):
    '''
    This function implements third-pass filter on extracted signals
    using their intrinsic LD relationships.
    
    Returns:
    --------
    clusters - list of Cluster objects
    '''
    import variants
    from utils import efficient_query
    import logging
    import pandas as pd
    import copy
    
    logger = logging.getLogger("root")
    
    ld_df = pd.read_csv(ld_fp, sep="\s+")
    contra_cluster_cnt = 0
    for cluster in clusters:

        if cluster.minp >= p_thresh:
            continue

        logger.debug(f"Current cluster: {cluster}")
        logger.debug(f"Significant variants in current cluster: {cluster.get_significant_vars(p_thresh)}")

        for var in cluster.get_significant_vars(p_thresh):

            new_cluster = copy.deepcopy(cluster)
            new_cluster.clear_r2()

            # Add omitted variants in extracted cluster coordinates
            add = efficient_query(all_snp_df, cluster.chrom, var.bp - 50000, var.bp + 50000)
            for _, _var in add.iterrows():
                new_cluster.add_variant(variants.Variant(_var.BP, _var.P, _var.CHR, _var.SNP))

            # Incorporate LD info
            new_cluster.add_r2_info(ld_df, var)

            # Evaluate p-value compatibility
            points = []
            var_cnt, contra_cnt = 0, 0
            for _var in new_cluster.variants:
                points.append(tuple([_var.r2, _var.logP]))
                if _var.r2 > 0.8 and _var.logP > 2:
                    var_cnt += 1
                elif _var.r2 > 0.8 and _var.logP <= 2:
                    contra_cnt += 1
            logger.debug(f"Number of variants with r2 > 0.8 and logP > 2: {var_cnt}")
            logger.debug(f"Contradictory evidence: {contra_cnt}")
            if var_cnt <= 1 or contra_cnt / (contra_cnt + var_cnt) > 0.8:
                clusters.remove(cluster)
                contra_cluster_cnt += 1
                break

    logger.info(f"Removed {contra_cluster_cnt} clusters with incompatibile p-values. ")
    
    return clusters

def classify_signals(clusters, p_thresh):
    
    from variants import Cluster
    import logging
    logger = logging.getLogger("root")
    
    singletons, duplets, real_clusters = [], [], []
    for cluster in clusters:
        n_signal = len(cluster.get_significant_vars(p_thresh))
        match n_signal:
            case 1:
                singletons.append(cluster)
            case 2:
                duplets.append(cluster)
            case 0:
                pass
            case default:
                real_clusters.append(cluster)
    logger.info(f"Extracted {len(singletons)} singletons, {len(duplets)} duplets and {len(real_clusters)} clusters. ")

    return singletons, duplets, real_clusters

def write_report(singletons, duplets, real_clusters, path, p_thresh, no_singleton=True):

    import os
    import pandas as pd
    import logging
    
    logger = logging.getLogger("root")
    
    report_out = []
    
    for cluster in real_clusters:
        report_out.append([cluster.chrom, cluster.start, cluster.end, (cluster.minp, cluster.maxp), len(cluster.get_significant_vars(p_thresh)), len(cluster.variants), "CLUSTER"])
    for duplet in duplets:
        report_out.append([duplet.chrom, duplet.start, duplet.end, (duplet.minp, duplet.maxp), len(duplet.get_significant_vars(p_thresh)), len(duplet.variants), "DUPLET"])
    if not no_singleton:
        for singleton in singletons:
            report_out.append([singleton.chrom, singleton.start, singleton.end, (singleton.minp, singleton.maxp), len(singleton.get_significant_vars(p_thresh)), len(singleton.variants), "SINGLETON"])
    report = pd.DataFrame(report_out)
    report.columns = ["chr", "start", "end", "p", "n_sig_variant", "n_variant", "type"]
    
    report.to_csv(os.path.join(path, "classification_report.tsv"), sep="\t", index=False)
    logger.info(f"Report written to {os.path.join(path, 'classification_report.tsv')}")
    
    return