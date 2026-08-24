# HiCBZIP

HiCBZIP is a Bayesian zero-inflated imputation framework for sparse single-cell Hi-C contact maps. This repository contains the method implementation and manuscript reproducibility code for:

1. The biologically realistic simulation study.
2. The downsampled NPC chrX recovery benchmark.
3. The SCORE mouse oocyte-to-zygote downstream embedding benchmark.

The repository is organized to help readers reproduce the main analyses from processed inputs, while keeping the focus on HiCBZIP code and final manuscript workflows. Large intermediate files and generated figures are not intended to be committed directly to GitHub.

## Repository Layout

- `HiCBZIP/`: core HiCBZIP model code, including the empirical-Bayes implementation, Stan models, and imputation functions.
- `simulation/`: scripts for the manuscript simulation study.
- `real_data_NPC_chrX/`: scripts for the downsampled NPC chrX benchmark.
- `real_data_SCORE_oocyte_zygote/`: scripts for the SCORE oocyte-to-zygote benchmark.
- `benchmarks/`: version and command notes for external comparison methods.
- `data/`: small example data and documentation for required processed inputs.
- `environment/`: software and package version information.
- `results/`: regenerated outputs; most files here should stay out of Git.
- `WORKFLOW.md`: ordered, study-by-study reproducibility map.

## Reproducibility Scope

This release is designed to reproduce manuscript-level results from processed inputs. It does not preserve every exploratory notebook or local benchmark runner used during method development. External methods such as scHiCluster, HiCImpute, Higashi, Fast-Higashi, and SCORE are documented by version/command notes and consumed as processed benchmark outputs where appropriate.

## Main Entry Points

- `examples/run_toy_all_HiCBZIP_variants.R`
- `simulation/03_run_HiCBZIP_GB_NB_simulation.R`
- `simulation/04_run_HiCBZIP_NGS_simulation.R`
- `simulation/05_run_HiCBZIP_NM_simulation.R`
- `simulation/06_prepare_simulation_manuscript_summaries.R`
- `real_data_NPC_chrX/01_make_NPC_chrX_processed_input.R`
- `real_data_NPC_chrX/02_run_HiCBZIP_GB_NB_NPC_chrX_all_coverage.R`
- `real_data_NPC_chrX/03_run_HiCBZIP_NM_NPC_chrX_one_coverage.R`
- `real_data_NPC_chrX/04_summarize_NPC_chrX_manuscript_metrics.R`
- `real_data_SCORE_oocyte_zygote/01_build_SCORE_HiCBZIP_inputs.R`
- `real_data_SCORE_oocyte_zygote/02_run_SCORE_embeddings.R`
- `real_data_SCORE_oocyte_zygote/03_summarize_SCORE_manuscript_metrics.R`

## Required Inputs

See `data/README.md` for the expected processed input layout.

## Public Source Data

The simulation regeneration scripts start from human brain snm3C-seq indexed-contact files from GSE130711, samples GSM3749700, GSM3750251, and GSM3751478. The two real-data studies use GSE72697/GSM1868576 for the NPC chrX recovery benchmark and GSE305523 plus GSE80006 for the SCORE oocyte-to-zygote benchmark.

## License

This project is released under the MIT License.

## Toy Example

The repository includes a small toy dataset:

```text
data/toy_chr5_diag20_cells8_lambda0.2.rds
```

To run the toy example from `code_clean`:

```r
source("examples/run_toy_all_HiCBZIP_variants.R")
```

This runs HiCBZIP-GB(NB), HiCBZIP-N(GS), and HiCBZIP-N(M), then writes one
output file and one 8-cell heatmap per method to `examples/`.

The two Stan variants use short toy MCMC settings by default; increase the
warmup and sampling iterations for real analyses.
