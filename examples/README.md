# Examples

Run the toy example from the `code_clean` root:

```r
source("examples/run_toy_all_HiCBZIP_variants.R")
```

The script loads `data/toy_chr5_diag20_cells8_lambda0.2.rds`, runs all three
HiCBZIP variants, and saves outputs in this folder:

```text
toy_raw_observed_cells8_heatmaps.png
toy_HiCBZIP_GB_NB_cells8_heatmaps.png
toy_HiCBZIP_NGS_cells8_heatmaps.png
toy_HiCBZIP_NM_cells8_heatmaps.png
toy_HiCBZIP_GB_NB_output.rds
toy_HiCBZIP_NGS_output.rds
toy_HiCBZIP_NM_output.rds
```

The Stan-based variants use toy MCMC defaults and silence Stan sampling output
through the `silent_stan` option in the HiCBZIP wrapper functions.
