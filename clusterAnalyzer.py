#!/usr/bin/env python3

import argparse
import os
import logging

import pandas as pd

from utils import *
from prepare_clusters import *

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--assoc', required = True, help='summary statistics from PLINK (*.logistic)')
    parser.add_argument('--bfile', required = False, help='PLINK binary file prefix')
    parser.add_argument('--min_rows', default = 10, help = 'Minimum SNP# requirement for a signal to be assessed')
    parser.add_argument('--mae', default = 1.0, help = 'Mean absolute error threshold for determining LD-p value consistency')
    parser.add_argument('--outdir', required = False, help='output directory (default: directory of summary statistics')
    parser.add_argument("--verbose", default = "INFO", help = "level of verbosity (default: INFO)")
    args = parser.parse_args()

    # Initialize a logger object
    logging.basicConfig(format='[%(asctime)s] %(levelname)s: %(message)s', datefmt='%a %b-%m %I:%M:%S%P', level = args.verbose.upper())
    logger = logging.getLogger("root")

    outd = args.outdir or os.path.dirname(args.assoc)
    os.makedirs(outd, exist_ok=True)

    # Extract significant SNPs
    assoc_df = pd.read_csv(args.assoc, sep="\s+")
    if "P" not in assoc_df.columns:
        raise ValueError("P is not found in columns. Are you using PLINK output?")
    selected_variants = assoc_df[assoc_df["P"].lt(0.05)].reset_index(drop=True)
    logger.info(f"Number of SNPs with p < 0.05: {selected_variants.shape[0]}")

    # Step 2: Second pass - Extract significant SNP clusters
    clusters, snps = [], []
    for chrom, sub_group in selected_variants.groupby("CHR"):
        filtered_dataframes = split_dataframe_by_cluster(sub_group.reset_index(drop=True))
        clusters.extend(filtered_dataframes)
        snps.extend([snp for cluster in filtered_dataframes for snp in cluster["SNP"]])

    # LD annotation
    ld_data = get_ld(outd, snps, args.bfile)
    ld_selected_clusters = filter_clusters(clusters, ld_data, assoc_df, outd, args.mae)

    # Reporting
    reporting(outd, assoc_df)
    logger.info("Cluster analysis completed successfully!")
