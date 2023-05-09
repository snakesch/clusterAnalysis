#!/usr/bin/env python3

def efficient_query(assoc_df, chrom, start, end):

    import pandas as pd

    search_start = assoc_df["CHR"].searchsorted(int(chrom))
    search_end = assoc_df["CHR"].searchsorted(int(chrom) + 1)

    by_chrom = assoc_df.iloc[search_start: search_end, :]
    pos_start = by_chrom.BP.searchsorted(start) + by_chrom.index[0]
    pos_end = by_chrom.BP.searchsorted(end) + by_chrom.index[0]

    return assoc_df.iloc[pos_start:pos_end+1, :]

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
        sig_var = len(cluster.get_significant_vars(p_thresh))
        tol_var = len(cluster.variants)
        report_out.append([cluster.chrom, cluster.start, cluster.end, (cluster.minp, cluster.maxp), sig_var, sig_var/tol_var])

    report = pd.DataFrame(report_out)
    report.columns = ["chr", "start", "end", "p", "n_signal", "frac_significant"]
    report = report.drop_duplicates(subset = ["chr", "start", "end"], keep = "first")

    report.to_csv(os.path.join(path, "cluster_summary.tsv"), sep="\t", index=False)
    logger.info(f"Cluster summary written to {os.path.join(path, 'cluster_summary.tsv')}")

    return
