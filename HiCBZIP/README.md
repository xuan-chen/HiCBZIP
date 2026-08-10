# HiCBZIP method code

Core files:

- `BHZIP_EB.R`: empirical Bayes Gamma-Beta HiCBZIP implementation, including neighborhood borrowing for HiCBZIP-GB(NB).
- `BHZIP_match_normal.stan`: Stan model for matched Normal-prior fitting used by HiCBZIP-N(M)-style analyses.
- `analysis_helpers.R`: matrix conversion, simulation, heatmap, SCC, and insulation-score helper functions used by manuscript workflows.

## Minimal example

Run the small example from the repository root:

```r
source("code_clean/HiCBZIP/BHZIP_EB.R")
source("code_clean/examples/run_minimal_HiCBZIP_example.R")
```

The example uses a tiny synthetic count matrix and demonstrates `run_BZIP_GB_NB()`.

## Dependencies

Core empirical Bayes functions use base R plus `dplyr`/`magrittr` for helper code. Matched Normal-prior fitting requires `cmdstanr` and a working CmdStan installation.

