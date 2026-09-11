# Benchmark Methods

This directory gives a top-level provenance summary for benchmark methods used in the HiCBZIP manuscript. Study-specific command workflows are stored with the corresponding study folder:

- `simulation/benchmark_methods/`
- `real_data_NPC_chrX/benchmark_methods/`
- `real_data_SCORE_oocyte_zygote/`

The main repository workflow focuses on HiCBZIP method code, archived manuscript-scale data objects, and final summary scripts. Full benchmark-method working directories, logs, temporary files, and large generated outputs are excluded from Git.

## Provenance Summary

| Method | Software/environment record | Used in |
| --- | --- | --- |
| HiCImpute | R package version 1.0. | Simulation study and NPC chrX recovery study. |
| scHiCluster | Python package `schicluster` 1.3.5.dev22+gd566046 was recorded in the `schicluster` conda environment. The `scHiCluster`/`schicluster` executable was not present on PATH in the recorded environment check; study scripts call the `hicluster impute-cell` interface used in the original workflow. | Simulation study, NPC chrX recovery study, and SCORE integrated-result comparison. |
| Higashi | Higashi 0.1.0a0 in the recorded `server_higashi_env_simulation` environment. | Simulation study, NPC chrX recovery study, and SCORE integrated-result comparison. |
| Fast-Higashi | Fast-Higashi 0.1.1a0 in the recorded `server_fast_higashi_env_simulation` environment. | Simulation study, NPC chrX recovery study, and SCORE integrated-result comparison. |
| SCORE | Python environment records are provided for SCORE and SCORE/Higashi-compatible runs. | SCORE mouse oocyte-to-zygote downstream embedding study. |

## Study-Specific Outputs

| Study | Archived or expected processed object | Notes |
| --- | --- | --- |
| Simulation | `data/processed/simulation/muS_combined_plot_260124.RData` | Combined matrix object used for simulation SMSE, SCC, and insulation-score summaries. |
| Simulation | `data/processed/simulation/higashi_unified_allchr_260208.RData` | Processed Higashi simulation output; the manuscript comparison uses the `nbr5` output. |
| Simulation | `data/processed/simulation/muS_combined_plot_both_fasthigashi_260613.RData` | Processed Fast-Higashi heatmap/summary input. |
| NPC chrX recovery | `data/processed/NPC_chrX/list_muS_unified_260614.RData` | Combined matrix object for Raw, HiCBZIP, HiCImpute, scHiCluster, Higashi, and Fast-Higashi across the reported downsampling levels. |
| SCORE oocyte-to-zygote | `data/processed/SCORE_oocyte_zygote/` and generated SCORE metric JSON folders | Inputs and integrated-result summaries used by the SCORE downstream benchmark. |

## Scope

Benchmark-method scripts in this repository are included to document the commands and file transformations used to produce the processed objects consumed by the manuscript summaries. They are not intended to replace the benchmark packages' own installation instructions.
