#!/usr/bin/env python3

## utils.py
## This script contains utility functions of ClusterAnalyzer.

import pandas as pd
import os
import logging

logger = logging.getLogger("root")

def get_ld(outd, snps, bfile):
    '''Writes LD information of SNPs to outd using genotypes stored in plink bfile'''
    import subprocess
    logger = logging.getLogger("root")
    
    if os.path.exists(ld_path := os.path.join(outd, os.path.basename(bfile).split('.')[0]) + ".ld"):
        logger.warning(f"Using previous LD results - {ld_path}")
        return pd.read_csv(ld_path, sep="\s+")
    
    # Write SNP list
    snp_list_path = os.path.join(outd, "snps.id")
    with open(snp_list_path, "w") as f:
        for snp in snps:
            f.write(snp)
            f.write("\n")
            
    subprocess.run(f"plink --bfile {bfile} --ld-snp-list {snp_list_path} --ld-window 10 --ld-window-r2 0.0 --r2 --out {os.path.join(outd, os.path.basename(bfile).split('.')[0])} 2>/dev/null", shell=True)
    
    ld_path = os.path.join(outd, os.path.basename(bfile).split('.')[0]) + ".ld"
    
    ld_data = pd.read_csv(ld_path, sep="\s+")
    logger.info(f"Loaded LD from {ld_path}")
    
    return ld_data

def reporting(outd, assoc_df):
    '''Creates summaries of prioritized clusters and variants based on association signals in assoc_df'''
    
    logger.info("Writing summary ... ")
    
    # Cluster reporting
    outp = os.path.join(outd, "merged_clusters.bed")
    cluster_df = pd.read_csv(outp, sep="\t", header=None)
    new_df = []
    var_cnt = 0
    variant_summary = open(os.path.join(outd, "variant_summary.assoc"), "w")
    headers = assoc_df.columns.tolist()
    for h in headers:
        variant_summary.write(h)
        variant_summary.write("\t")
    variant_summary.write("\n")
    for _, region in cluster_df.iterrows():
        chrom, start, end = region
        assoc_by_chrom = assoc_df[assoc_df["CHR"] == chrom]
        assoc_target = assoc_by_chrom[assoc_by_chrom["BP"].between(start, end)]
        var_cnt += assoc_target.shape[0]
        assoc_target.to_csv(os.path.join(outd, "variant_summary.assoc"), sep="\t", index=False, header=False, mode="a")
        region["SignalCount"] = "{}/{}".format(assoc_target[assoc_target["OR"] > 1].shape[0], assoc_target[assoc_target["OR"] < 1].shape[0])
        region["Pmin"] = "{}/{}".format(assoc_target[assoc_target["OR"] > 1]["P"].min(), assoc_target[assoc_target["OR"] < 1]["P"].min())
        region["Pmax"] = "{}/{}".format(assoc_target[assoc_target["OR"] > 1]["P"].max(), assoc_target[assoc_target["OR"] < 1]["P"].max())
        new_df.append(region)

    reportCluster = pd.DataFrame(new_df)
    logger.info(f"Number of prioritized clusters: {reportCluster.shape[0]}")
    logger.info(f"Number of variants in clusters: {var_cnt}")
    reportCluster.columns = headers
    reportCluster.to_csv(os.path.join(outd, "cluster_summary.tsv"), sep="\t", header=True, index=False)
                
    return
                
def cleanup(outd, bfile):
    '''Cleans up intermediate files in outd'''
    
    snp_list = os.path.join(outd, "snps.id")    
    ld_path = os.path.join(outd, os.path.basename(bfile) + ".ld")
    merged_bed = os.path.join(outd, "merged_clusters.bed")
                
    for file in [snp_list, ld_path, merged_bed]:
        if os.path.exists(file):
            os.remove(file)
                
    return