#!/usr/bin/env python3

## prepare_clusters.py
## This script contains functions used to determine genomic clusters based on association signals of summary statistics; min_rows determines how many SNPs required for a cluster to be considered to have sufficient points for analysis.

import pandas as pd
import numpy as np
import os

def split_dataframe_by_cluster(df, threshold=5000, min_rows=10):
    '''Splits a DataFrame into clusters by genomic distances and ORs'''
    
    if 'BP' not in df.columns:
        raise ValueError("Column 'BP' does not exist in the DataFrame.")
    
    # Sort the DataFrame by the 'BP' column
    df_sorted = df.sort_values('BP')
    
    # Identify the split indices
    split_indices = []
    prev_value = None
    for index, row in df_sorted.iterrows():
        value = row['BP']
        if prev_value is not None and value - prev_value > threshold:
            split_indices.append(index)
        prev_value = value
    
    # Split the DataFrame based on the identified indices
    split_dataframes = []
    start_index = 0
    for split_index in split_indices:
        split_df = df_sorted.iloc[start_index:split_index]
        if len(split_df) > min_rows:
            split_dataframes.append(split_df)
        start_index = split_index
        
    return split_dataframes

def filter_clusters(clusters, ld_data, assoc_df, outd, mae):
    '''Determine which clusters have consistent p-values with LD.'''
    
    from sklearn.metrics import mean_absolute_error
    from pybedtools import BedTool as bt
    
    selected_clusters = []

    for cluster in clusters:
        start, end = cluster["BP"].tolist()[0], cluster["BP"].tolist()[-1]
        chrom = cluster["CHR"].tolist()[0]
        
        # Select the most significant SNP
        most_significant_snp = cluster.sort_values("P").iloc[0, :]
        
        # Extract LD of related SNPs
        ld_target = ld_data[(ld_data["CHR_A"] == chrom) & (ld_data["BP_A"] == most_significant_snp["BP"])]
        linked_snps = ld_target[ld_target["R2"] > 0.4]
        
        # Extract OR and p values from summary statistics file
        cluster_dict_or_p = { b: (o, p) for b, o, p in zip(cluster["BP"].tolist(), cluster["OR"].tolist(), cluster["P"].tolist()) }
        
        # Identify good points based on LD
        good_points = {}
        for bp, assoc in cluster_dict_or_p.items():
            target = linked_snps[linked_snps["BP_B"] == bp]
            if target.shape[0] == 0:
                continue
            good_points[bp] = (target["R2"].tolist()[0], assoc[0], assoc[1])
        
        good_points_df = pd.DataFrame.from_dict(good_points, orient="index")
        good_points_df.columns = ["R2", "OR", "P"]
        good_points_df["-logP"] = -np.log(good_points_df["P"])
        signals = good_points_df[good_points_df["OR"] > 1], good_points_df[good_points_df["OR"] < 1]
        
        for signal in signals:
            if signal.shape[0] < 5:
                continue
            else:
                beta = signal.sort_values("-logP")["-logP"].tolist()[-1]
                signal.loc[:, "Predicted association"] = beta * signal["R2"]
                
                # Evaluate MAE
                mar = mean_absolute_error(signal["-logP"], signal["Predicted association"])

                # Keep cluster if MAE <= 3.0
                if mar <= 3.0:
                    selected_clusters.append((chrom, start, end))
                    break # No need to evaluate the other signal if already prioritized
    
    # Merge overlapping clusters
    df = pd.DataFrame(selected_clusters)
    c = bt.from_dataframe(df)
    outp = os.path.join(outd, "merged_clusters.bed")
    c.sort().merge().saveas(outp)
    
    return pd.read_csv(outp, sep="\t", header=None)