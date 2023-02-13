#!/usr/bin/env python3

import argparse
import os
from datetime import datetime as dt
import logging
import pandas as pd
import variants
import utils
import subprocess

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--snp_df', required = True, help='summary statistics from PLINK (*.logistic)')
    parser.add_argument('--p_thresh', required = True, type = float, help='p-value threshold for significance classification')
    parser.add_argument('--bfile', required = True, help='PLINK binary file prefix (with full path)')
    parser.add_argument("--no_singletons", action = "store_true", help = "do not report singletons")
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
    window = args.window
    report_all = args.no_singletons
    
    wkd = os.path.dirname(assoc_fp)
    ld_fp = os.path.join(wkd, "run_" + dt.today().strftime("%Y%m%d") + ".ld")
    
    # Load summary statistics
    assoc_df = pd.read_csv(assoc_fp, sep="\s+")
    if "P" not in assoc_df.columns:
        raise ValueError("P is not found in columns. Are you using PLINK output?")
    
    # Step 1: First pass - Filter out SNP with p >= 0.05
    assoc_sig = assoc_df[assoc_df["P"].lt(0.05)]
    logger.info(f"Number of signals with p < 0.05 : {assoc_sig.shape[0]}")

    # Step 2: Second pass - Filter out clusters with only one outlier & generate a SNP list of variants with p < p_thresh
    clusters = utils.second_pass_singletons(assoc_sig, wkd, window=window, p_thresh=p_thresh)

    # Step 3: Third pass - Use LD information to filter out signals with incompatible p-values
    subprocess.run(["plink", "--bfile", bfile_prefix, "--ld-snp-list", os.path.join(wkd, "snps.id"), "--ld-window", "99999", "--ld-window-kb", "100",
                   "--ld-window-r2", "0.0", "--out", ld_fp[:-3], "--r2"])
    clusters = utils.third_pass_ld(ld_fp, assoc_df, clusters, p_thresh)

    # Step 4: Signal classification and reporting
    singletons, duplets, real_clusters = utils.classify_signals(clusters, p_thresh)
    utils.write_report(singletons, duplets, real_clusters, wkd, no_singleton=report_all)
    
    logger.info("Cluster analysis completed successfully!")