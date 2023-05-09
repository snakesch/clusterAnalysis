# ClusterAnalyzer
ClusterAnalyzer is a python implementation of a novel clustering approach that filters out potential false positives due to genotyping error. ClusterAnalyzer generates a list of clusters which more likely represent loci of true association. As of its current implementation, it is intended for association studies with small power where substantial amounts of false positives may exist. Users are therefore reminded to interpret its result with caution. 

## Prerequisite
* Python v3.9
* NumPy v1.23.4
* Pandas v1.5.1

## Usage

```{bash}
usage: clusterAnalyzer.py [-h] --snp_df SNP_DF [--bfile BFILE] [--ld LD] [--p_thresh P_THRESH] [--outdir OUTDIR] [--verbose VERBOSE]

optional arguments:
  -h, --help           show this help message and exit
  --snp_df SNP_DF      summary statistics from PLINK (*.logistic)
  --bfile BFILE        PLINK binary file prefix (with full path)
  --ld LD              pre-computed PLINK LD file
  --p_thresh P_THRESH  p-value threshold for SNP reporting
  --outdir OUTDIR      output directory (default: directory of summary statistics
  --verbose VERBOSE    level of verbosity (default: INFO)
```

## Result
ClusterAnalyzer reports extracted clusters in two files: `filtered_clusters.assoc` and `cluster_summary.tsv`. `filtered_clusters.assoc` includes summary statistics of variants found within reported genomic clusters in PLINK format; `cluster_summary.tsv` includes various metrics of reported clusters which are explained further below.

| chr  | start | end | p | n_signal | frac_significant |
| ---- | ------------- | ------------- | ------------- | -------------- | ------------ |
| 1  | 163188121 | 163200062 | (0.000159, 0.03013) | 7 | 0.0345 |
| 1  | 163241577 | 163284768 | (0.0001622, 0.03494) | 8 | 0.0547 |

* `chr`: the chromosome where the signal is found
* `start`: start coordinate of the signal
* `end`: end coordinate of the signal
* `p`: p-value range (min p, max p) of variants in the cluster
* `n_signal`: number of variants with p < p_thresh in the cluster
* `frac_significant`: fraction of significant SNPs in the cluster

## Benchmark


## Correspondence
Louis SHE (snakesch@connect.hku.hk)

