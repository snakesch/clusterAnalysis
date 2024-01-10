# ClusterAnalyzer
ClusterAnalyzer is a clustering algorithm for prioritizing likely true signals from an excess of false positives. The tool is intended for genome-wide association studies with minimal sample size which suffer from insufficient power. Users are reminded to interpret the results with caution.

## Prerequisite
* Python >= v3.9
* NumPy >= v1.23.4
* Pandas >= v1.5.1

## Usage
```{bash}
usage: clusterAnalyzer.py [-h] --assoc ASSOC [--bfile BFILE] [--min_rows MIN_ROWS] [--mae MAE] [--outdir OUTDIR]
                          [--verbose VERBOSE]

options:
  -h, --help           show this help message and exit
  --assoc ASSOC        summary statistics from PLINK (*.logistic)
  --bfile BFILE        PLINK binary file prefix
  --min_rows MIN_ROWS  Minimum SNP# requirement for a signal to be assessed
  --mae MAE            Mean absolute error threshold for determining LD-p value consistency
  --outdir OUTDIR      output directory (default: directory of summary statistics
  --verbose VERBOSE    level of verbosity (default: INFO)
```

## Input/Output
ClusterAnalyzer takes association summary statistics as input and prioritizes signals that are likely to be true. Prioritized results are written to OUTDIR with variant_summary.assoc containing individual SNPs and cluster_sumamry.txt containing all prioritized genomic regions.  

| Chr  | Start | End | SignalCount | Pmin | Pmax |
| ---- | ------------- | ------------- | ------------- | -------------- | ------------ |
| 1  | 163188121 | 163200062 | 15/9 | 0.0003528/0.03249 | 0.9963/0.9837 |
| 1  | 163241577 | 163284768 | 78/99 | 0.0005863/0.03234 | 0.9819/0.9716 |

* `Chr`: the chromosome where the signal is found
* `Start`, `End`: start and end coordinates of the signal
* `SignalConut`: number of SNPs with OR > 1 and number of SNPs with OR < 1 in the cluster
* `Pmin`, `Pmax`: minimum and maximum p values of SNPs with OR > 1 and OR < 1 in each cluster

## Code correspondence
Louis SHE (snakesch@connect.hku.hk)

