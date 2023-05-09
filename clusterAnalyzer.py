#!/usr/bin/env python3

import argparse
import os
import subprocess
import logging
import copy
from datetime import datetime as dt

import numpy as np
import pandas as pd

from variants import Variant, Cluster
from utils import *

def extract_analytic_clusters(assoc_sig, prom_thresh, fout, window=5000):

    '''
    Extracts plausible genomic loci based on number of outliers

    assoc_sig: dataframe of SNP p < 0.05
    outlier_thresh: outlier threshold for selecting genomc loci
    fout: output fs

    '''
    logger = logging.getLogger("root")

    def extract_clusters(assoc_grp_df, window):

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

def assess_correlation(clusters, assoc_df, ld_df, resid_cutoff=0.5, prom_thresh=5, ld_cutoff=0.4, plot=False):

    '''
    Inputs a complete of candidate clusters, a df of all SNPs, coef_cutoff: cutoff of
    correlation coef (R2; default = 0.5), prom_thresh: number of similarly significant
    SNPs required. ld_cutoff: LD r2 cutoff above which degree of fitness to y = ymax * x
    is assessed.

    Returns:
    --------
    filtered_clusters

    '''
    def getMAR(points, min_size = 3):

        '''
        This function takes a numpy array (r2, logP) and computes the respective mean absolute residual (MAR).

        Return:
        --------
        float: mean absolute residual
        999: if insufficient points

        '''

        if points.shape[0] < min_size or points[points[:, 0] == 1.0].shape[0] == 0:
            return 999

        to_fit = points[points[:, 0] > 0.4]
        refp = points[points[:, 0] == 1.0, 1][0]

        X = to_fit[:, 0].reshape(-1, 1)
        y_true = to_fit[:, 1]
        y_pred = refp * X
 
        return np.mean(np.absolute(y_true - y_pred))

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
        new_cluster.add_r2_info(ld_df, cluster.variants[0])
        points = np.array([ (_var.r2, _var.logP, _var.bp) for _var in new_cluster.variants ])

        # Keep track of processed SNPs (Exclude tightly linked SNPs)
        linked_var = np.empty(1)
        while len(new_cluster.variants) > prom_thresh:

            # Initialize the points
            for var in new_cluster.variants:
                if var.logP == points[points[:, 1].argsort()[::-1]][0, 1]:
                    break
            new_cluster.add_r2_info(ld_df, var)
            points = np.array([ (_var.r2, _var.logP, _var.bp) for _var in new_cluster.variants ])

            # Case: No LD information available
            if points[points[:, 0] == 1.0].shape[0] == 0:
                linked_var = np.append(linked_var, points[0, 2])
            # Case: Insufficient points for fitting
            elif points[points[:, 0] > 0.8].shape[0] < prom_thresh:
                linked_var = np.append(linked_var, points[points[:, 0] > 0.8, 2])
            else:
                # Compute MAR
                linked_var = np.append(linked_var, points[points[:, 0] > 0.8, 2])
                if getMAR(points, min_size = prom_thresh) < resid_cutoff:
                    filtered_clusters.append(copy.deepcopy(new_cluster))
                    break

            _linked = list(map(int, linked_var))
            new_cluster.variants = sorted([var for var in new_cluster.variants if var.bp not in _linked],
                                          key = lambda var: var.logP, reverse = True)
            points = points[points[:, 0] <= 0.8, :]
            new_cluster.clear_r2()
            if points.shape[0] < prom_thresh:
                break

    return filtered_clusters

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--snp_df', required = True, help='summary statistics from PLINK (*.logistic)')
    parser.add_argument('--bfile', required = False, help='PLINK binary file prefix (with full path)')
    parser.add_argument('--ld', required = False, help='pre-computed PLINK LD file')
    parser.add_argument('--p_thresh', default = 0.05, type = float, help='p-value threshold for SNP reporting')
    parser.add_argument('--window', default = 5000, help='window size for finding signal aggregates (default: 5000)')
    parser.add_argument('--outdir', required = False, help='output directory (default: directory of summary statistics')
    parser.add_argument("--verbose", default = "INFO", help = "level of verbosity (default: INFO)")
    args = parser.parse_args()

    # Initialize a logger object
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P',
                            level = args.verbose.upper())
    logger = logging.getLogger("root")

    assoc_fp = args.snp_df
    p_thresh = args.p_thresh
    bfile_prefix = args.bfile

    PROM_THRESH = 5
    RESID_CUTOFF = 1.0
    LD_CUTOFF = 0.4

    wkd = args.outdir or os.path.dirname(assoc_fp)
    ld_fp = args.ld or os.path.join(wkd, "run_" + dt.today().strftime("%Y%m%d") + ".ld")

    # Load summary statistics
    assoc_df = pd.read_csv(assoc_fp, sep="\s+")
    if "P" not in assoc_df.columns:
        raise ValueError("P is not found in columns. Are you using PLINK output?")

    # Step 1: First pass - Filter out SNP with p >= 0.05
    assoc_sig = assoc_df[assoc_df["P"].lt(0.05)].reset_index(drop=True)
    logger.info(f"Number of signals with p < 0.05 : {assoc_sig.shape[0]}")

    # Step 2: Second pass - Extract significant SNP clusters
    fout = open(os.path.join(wkd, "snps.id"), "w")
    clusters = extract_analytic_clusters(assoc_sig, PROM_THRESH, fout, window=args.window)
    logger.info(f"SNP ID list written to {os.path.join(wkd, 'snps.id')} ")
    fout.close()

    # Step 3: Third pass - Extract SNP clusters with consistent LD
    if not args.ld:
        if not args.bfile:
            raise FileNotFoundError("No LD data provided. Stop. ")
        else:
            subprocess.run(["plink", "--bfile", bfile_prefix, "--ld-snp-list", os.path.join(wkd, "snps.id"), "--ld-window", "99999", "--ld-window-kb", "100",
                           "--ld-window-r2", "0.0", "--out", ld_fp[:-3], "--r2"])

    logger.info("Loading LD file ... ")
    ld_df = pd.read_csv(ld_fp, sep="\s+")
    logger.info("Loaded LD file.")
    logger.info("Begin to evaluate LD consistency ... ")
    filtered_clusters = assess_correlation(clusters, assoc_df, ld_df, resid_cutoff = RESID_CUTOFF, prom_thresh = PROM_THRESH, ld_cutoff = LD_CUTOFF)

    # Step 4: Reporting
    write_ss(assoc_df, filtered_clusters, os.path.join(wkd, "filtered_clusters.assoc"), p_thresh = p_thresh)
    write_summary(filtered_clusters, wkd, p_thresh = 0.05)

    logger.info("Cluster analysis completed successfully!")
