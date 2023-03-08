#!/usr/bin/env python3

def efficient_query(assoc_df, chrom, start, end):

    import pandas as pd

    search_start = assoc_df["CHR"].searchsorted(int(chrom))
    search_end = assoc_df["CHR"].searchsorted(int(chrom) + 1)

    by_chrom = assoc_df.iloc[search_start: search_end, :]
    pos_start = by_chrom.BP.searchsorted(start) + by_chrom.index[0]
    pos_end = by_chrom.BP.searchsorted(end) + by_chrom.index[0]

    return assoc_df.iloc[pos_start:pos_end+1, :]

def extract_analytic_clusters(assoc_sig, prom_thresh, fout, window=5000):

    '''
    Extracts plausible genomic loci based on number of outliers

    assoc_sig: dataframe of SNP p < 0.05
    outlier_thresh: outlier threshold for selecting genomc loci
    fout: output fs

    '''

    import logging
    import pandas as pd
    import numpy as np
    from variants import Variant, Cluster

    logger = logging.getLogger("root")

    def extract_clusters(assoc_grp_df, window):

        import pandas as pd
        import numpy as np

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

    clusters, singletons = [], 0
    sig_var_cnt = 0
    for chrom, grp_df in assoc_sig.groupby("CHR"):
        grp_df = grp_df.reset_index(drop=True)
        _clusters = extract_clusters(grp_df, window)
        for _cand in _clusters:
            var_lst = [  Variant(var.BP, var.P, var.CHR, var.SNP)
                          for idx, var in grp_df.iloc[list(_cand), :].iterrows()  ]
            cluster_cand = Cluster(var_lst)
            n_outlier = cluster_cand.get_outliers()

            if n_outlier <= prom_thresh and n_outlier != 0:
                continue
            elif n_outlier == 0:
                mean = sum([_var.logP for _var in cluster_cand.variants])/len(cluster_cand.variants)
                percentiles = np.percentile([_var.logP for _var in cluster_cand.variants], [25, 75])
                upper_bdd = mean + (percentiles[1] - percentiles[0]) * 1.5
                if upper_bdd <= 2.5:
                    continue
                else:
                    clusters.append(cluster_cand)
                    sig_var_cnt += len(var_lst)
                    for _var in var_lst:
                        fout.write(_var.snpid)
                        fout.write("\n")
            else:
                clusters.append(cluster_cand)
                sig_var_cnt += len(var_lst)
                for _var in var_lst:
                    fout.write(_var.snpid)
                    fout.write("\n")

    logger.info(f"Extracted {sig_var_cnt} SNPs in {len(clusters)} clusters. ")

    return clusters

def assess_correlation(clusters, assoc_df, ld_df, resid_cutoff=0.5, prom_thresh=5, ld_cutoff=0.4):

    '''
    Inputs a complete of candidate clusters, a df of all SNPs, coef_cutoff: cutoff of
    correlation coef (R2; default = 0.5), prom_thresh: number of similarly significant
    SNPs required. ld_cutoff: LD r2 cutoff above which degree of fitness to y = ymax * x
    is assessed.

    Returns:
    --------
    filtered_clusters

    '''

    import copy
    import numpy as np
    from variants import Variant, Cluster
    import logging

    def get_residuals(points, min_size = 3):

        '''
        This function takes a numpy array of points (r2, logP) and computes the mean residuals
        based on a linear model of y_hat = y_ref * r2, where y_ref is the logP value of LD re-
        ference point.

        Return:
        --------
        float: a non-negative mean residual estimate

        '''
        import numpy as np

        to_fit = points[points[:, 0] > 0.4]
        refp = points[points[:, 0] == 1.0, 1][0]

        # Too few signals to fit a linear model
        if to_fit.shape[0] < min_size:
            return 999

        X = to_fit[:, 0].reshape(-1, 1)
        y_true = to_fit[:, 1]

        y_pred = refp * X
        resid = np.mean(np.absolute(y_true - y_pred))

        return resid # if resid_var < 1.5 else resid_var

    logger = logging.getLogger("root")

    filtered_clusters = []
    for cluster in clusters:

        cluster.variants = sorted(cluster.variants, key = lambda x: x.logP, reverse=True)
        new_cluster = copy.deepcopy(cluster)

        # Add omitted variants in extracted cluster coordinates
        add = efficient_query(assoc_df, cluster.chrom, cluster.start - 50000, cluster.end + 50000)
        for _, _var in add.iterrows():
            if _var.P >= 0.05:
                new_cluster.add_variant(Variant(_var.BP, _var.P, _var.CHR, _var.SNP))

        # Keep track of processed SNPs (Exclude tightly linked SNPs)
        processed = []
        for var in cluster.variants:

            if var in processed:
                continue

            new_cluster.clear_r2()

            # Incorporate LD info
            new_cluster.add_r2_info(ld_df, var)

            # Evaluate p-value compatibility
            points = np.array([ (_var.r2, _var.logP) for _var in new_cluster.variants ])
            for _var in new_cluster.variants:
                if _var.r2 > 0.8:
                    processed.append(_var)

            resid = get_residuals(points, min_size = prom_thresh)

            # Case: Too few signals to fit linear model
            if resid == 999:
                new_cluster.variants.remove(var)
                continue
            elif resid < resid_cutoff and len(processed) >= prom_thresh:
                filtered_clusters.append(new_cluster) # with r2 info

    # Output the number of unique clusters
    regions = [ str(_cluster.chrom) + ":" + str(_cluster.start) for _cluster in filtered_clusters]
    logger.info(f"A total of {len(set(regions))} clusters passed all filters. ")

    return filtered_clusters

def plot_cluster(cluster, ld_cutoff=0.4):

    from matplotlib import pyplot as plt
    import numpy as np

    color_dict = {1.: "tab:purple"}

    points = np.array([ (_var.r2, _var.logP) for _var in cluster.variants ])

    fig, ax = plt.subplots()
    for x, y in zip(points[:, 0], points[:, 1]):
        ax.scatter(x, y, color = color_dict.get(x, "tab:blue"))
    ax.set_xlabel("LD $r^2$")
    ax.set_ylabel("$-logP$")
    ax.set_title(f"""Cluster {cluster.chrom}:{cluster.start}-{cluster.end}"""
                 f""" (LD ref: {cluster.ld_ref.chrom}:{cluster.ld_ref.bp})""")
    ax.axvline(x=ld_cutoff, c="0.4", linestyle = "--")
    refp = points[points[:, 0] == 1.0, 1][0]
    ax.plot(points[:, 0], points[:, 0] * refp, color = "tab:red")
    ax.set_xlim(left=-0.05, right=1.05)
    ax.set_ylim(bottom=-0.2, top=max(points[:, 1])+0.2)
    plt.show()

    return

def write_ss(assoc_df, clusters: list, outfp, p_thresh = 1.0):

    import pandas as pd
    import logging

    logger = logging.getLogger("root")

    clusters_coord = set()
    for cluster in clusters:
        clusters_coord.add((cluster.chrom, cluster.start, cluster.end))

    out_dfs = []
    for cluster in clusters_coord:
        out_dfs.append(efficient_query(assoc_df, cluster[0], cluster[1], cluster[2]))
    out_df = pd.concat(out_dfs, axis = 0)
    out_df = out_df[out_df["P"].le(p_thresh)]
    out_df = out_df.sort_values(by=["CHR", "BP"])
    out_df.to_csv(outfp, sep="\t", index = False)
    logger.info(f"Summary statistics written to {outfp} .")

    return out_df

def write_summary(clusters, path, p_thresh = 0.05):

    import os
    import pandas as pd
    import logging

    logger = logging.getLogger("root")

    report_out = []

    for cluster in clusters:
        report_out.append([cluster.chrom, cluster.start, cluster.end, (cluster.minp, cluster.maxp), len(cluster.get_significant_vars(p_thresh))])

    report = pd.DataFrame(report_out)
    report.columns = ["chr", "start", "end", "p", "n_signal"]

    report.to_csv(os.path.join(path, "cluster_summary.tsv"), sep="\t", index=False)
    logger.info(f"Cluster summary written to {os.path.join(path, 'summary.tsv')}")

    return
