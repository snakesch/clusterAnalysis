# ClusterAnalyzer
ClusterAnalyzer is a python implementation of a novel approach called cluster analysis described [here](). Essentially, ClusterAnalyzer attempts to exclude common variants that are found significant by genotyping error given small sample sizes; it prioritizes a list of clusters/duplets that are less likely to have been observed by genotyping error. As ClusterAnalyzer is intended for small samples, users are reminded that any results generated from ClusterAnalyzer shall only be interpreted as "suggestive association" and shall, by no means, be interpreted as of approximate power as standard GWAS analysis which employs a much larger cohort.

## Prerequisite
* Python v3.10 (or newer)
* NumPy v1.23.4
* Pandas v1.5.1

## Usage

```{bash}
usage: clusterAnalyzer.py [-h] --snp_df SNP_DF --p_thresh P_THRESH --bfile BFILE [--no_singletons] [--window WINDOW]
                          [--verbose VERBOSE]

options:
  -h, --help           show this help message and exit
  --snp_df SNP_DF      summary statistics from PLINK (*.logistic)
  --p_thresh P_THRESH  p-value threshold for significance classification
  --bfile BFILE        PLINK binary file prefix (with full path)
  --no_singletons      do not report singletons
  --window WINDOW      window size for finding signal aggregates (default: 5000)
  --verbose VERBOSE    level of verbosity (default: INFO)
```
Warning: Users are advised not to change `--window` as it may alter the expected behavior of ClusterAnalyzer.

## Result
ClusterAnalyzer extracts clusters, duplets that may possibly be associated with the phenotype of interest. Singletons can be excluded from the report with `--no_singletons`. An example of output report is provided below.

| chr  | start | end | p | n_sig_variant | n_variant | type |
| ---- | ------------- | ------------- | ------------- | ----------- | ------------- | ------------- |
| 1  | 163188121 | 163200062 | (0.000159, 0.03013) | 7 | 15 | CLUSTER |
| 1  | 163241577 | 163284768 | (0.0001622, 0.03494) | 8 | 32 | CLUSTER |

* `chr`: the chromosome where the signal is found
* `start`: start coordinate of the signal
* `end`: end coordinate of the signal
* `p`: p-value interval (min p, max p) of variants within the signal region
* `n_sig_variant`: number of variants with p < p_thresh within the signal region
* `n_variant`: number of variants with p < 0.05 within the signal region
* `type`: classified signal type (CLUSTER/DUPLET/SINGLETON)

## Correspondence
Louis SHE (snakesch@connect.hku.hk)
