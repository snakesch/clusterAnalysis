#!/usr/bin/env python3

import argparse
import os
import subprocess
import logging
import pandas as pd

from datetime import datetime as dt
from variants import Variant, Cluster
from utils import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--snp_df', required = True, help='summary statistics from PLINK (*.logistic)')
    parser.add_argument('--bfile', required = True, help='PLINK binary file prefix (with full path)')
    parser.add_argument('--p_thresh', default = 0.05, type = float, help='p-value threshold for SNP reporting')
    parser.add_argument('--window', default = 5000, help='window size for finding signal aggregates (default: 5000)')
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
    
    wkd = os.path.dirname(assoc_fp)
    ld_fp = os.path.join(wkd, "run_" + dt.today().strftime("%Y%m%d") + ".ld")
    
    # Load summary statistics
    assoc_df = pd.read_csv(assoc_fp, sep="\s+")
    if "P" not in assoc_df.columns:
        raise ValueError("P is not found in columns. Are you using PLINK output?")

    # Step 1: First pass - Filter out SNP with p >= 0.05
    assoc_sig = assoc_df[assoc_df["P"].lt(0.05)].reset_index(drop=True)
    logger.info(f"Number of signals with p < 0.05 : {assoc_sig.shape[0]}")

    # Step 2: Extract signal clusters and significant SNPs from within
    fout = open(os.path.join(wkd, "snps.id"), "w")
    clusters = extract_analytic_clusters(assoc_sig, PROM_THRESH, fout, window=args.window)
    logger.info(f"SNP ID list written to {os.path.join(wkd, 'snps.id')} ")
    fout.close()
    subprocess.run(["plink", "--bfile", bfile_prefix, "--ld-snp-list", os.path.join(wkd, "snps.id"), "--ld-window", "99999", "--ld-window-kb", "100",
                   "--ld-window-r2", "0.0", "--out", ld_fp[:-3], "--r2"])

    # Step 3: Compute LD r2 for reported SNPs
    logger.info("Loading LD file ... ")
    ld_df = pd.read_csv(ld_fp, sep="\s+")
    logger.info("Loaded LD file.")
    filtered_clusters = assess_correlation(clusters, assoc_df, ld_df, resid_cutoff = RESID_CUTOFF, prom_thresh = PROM_THRESH, ld_cutoff = LD_CUTOFF)

    # Step 4: Reporting
    write_ss(assoc_df, filtered_clusters, os.path.join(wkd, "filtered_clusters.assoc"), p_thresh = p_thresh)
    write_summary(filtered_clusters, wkd, p_thresh = 0.05)
    
    logger.info("Cluster analysis completed successfully!")