# HiCBZIP method code

Core files:

- `BHZIP_EB.R`: empirical Bayes Gamma-Beta HiCBZIP implementation, including neighborhood borrowing for HiCBZIP-GB(NB).
- `BHZIP_GM_nocov.stan`: Stan model for HiCBZIP-N(GS).
- `BHZIP_match_normal.stan`: Stan model for matched Normal-prior fitting used by HiCBZIP-N(M)-style analyses.
- `analysis_helpers.R`: matrix conversion, simulation, heatmap, SCC, and insulation-score helper functions used by manuscript workflows.

## Toy example

Run the small toy example from the `code_clean` root:

```r
source("examples/run_toy_all_HiCBZIP_variants.R")
```

The example loads `data/toy_chr5_diag20_cells8_lambda0.2.rds`, runs
HiCBZIP-GB(NB), HiCBZIP-N(GS), and HiCBZIP-N(M), and writes one output file and
one 8-cell heatmap per method to `examples/`.

## Dependencies

Core empirical Bayes functions use base R plus `dplyr`/`magrittr` for helper code. Matched Normal-prior fitting requires `cmdstanr` and a working CmdStan installation.

The Stan wrappers accept `silent_stan = TRUE` to suppress CmdStan sampling
output in lightweight examples or batch runs.
