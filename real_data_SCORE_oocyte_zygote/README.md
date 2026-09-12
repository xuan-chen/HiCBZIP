# Real-Data Study 2: SCORE Oocyte-To-Zygote Benchmark

This folder contains the public workflow for the mouse oocyte-to-zygote downstream benchmark. The analysis uses SCORE-formatted 1 Mb `mm10` input data and compares external-imputation pipelines against integrated embedding methods.

## Main Workflow

1. Place SCORE processed raw inputs under `data/processed/SCORE_oocyte_zygote/`.
2. Run `01_build_SCORE_HiCBZIP_inputs.R`.
   - Builds Raw, HiCBZIP-GB/GB(NB), and HiCBZIP-N(M) SCORE-ready `.scool` inputs.
3. Build or provide the scHiCImpute `.scool` input described in `benchmark_methods/README.md`.
4. Run `02_run_SCORE_embeddings.R`.
   - Runs SCORE InnerProduct and SnapATAC/no-IDF for Raw, HiCBZIP-GB/GB(NB), HiCBZIP-N(M), and scHiCImpute.
5. Run or provide the integrated-method SCORE outputs described in `benchmark_methods/README.md`.
   - scHiCluster, Higashi, and Fast-Higashi.
6. Run `03_summarize_SCORE_manuscript_metrics.R`.
   - Generates manuscript metric tables and figures from the SCORE result folders.

## Included Source Workflows

| Script | Purpose |
| --- | --- |
| `build_SCORE_HiCBZIP_GB_inputs.Rmd` | Builds Raw and HiCBZIP-GB/GB(NB) SCORE-ready pair files and `.scool` inputs. |
| `build_SCORE_HiCBZIP_NM_inputs.R` | Builds HiCBZIP-N(M) SCORE-ready pair files and `.scool` inputs. |
| `run_SCORE_innerproduct_four_inputs_10runs.R` | Runs SCORE InnerProduct over Raw, HiCBZIP-GB/GB(NB), HiCBZIP-N(M), and scHiCImpute. |
| `run_SCORE_snapatac_four_inputs_10runs.R` | Runs SCORE SnapATAC/no-IDF over the same four external-imputation inputs. |
| `summarize_SCORE_final_metrics.Rmd` | Final manuscript metric summary, PCA panel, and example heatmap workflow. |

Study-specific benchmark-method scripts are stored in `benchmark_methods/`.
