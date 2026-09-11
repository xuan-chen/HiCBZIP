# Real-data study 1: NPC chrX recovery

This workflow reproduces the downsampled chromosome-X recovery study using the NPC Hi-C matrix from GSE72697, sample GSM1868576.

Large processed inputs and generated outputs are not tracked in Git. Place them under the paths documented in `../data/README.md` before running the manuscript-scale workflow.

## Required Inputs

| File | Purpose |
| --- | --- |
| `data/processed/NPC_chrX/NPC250k_0h_X.mat` | Processed NPC chromosome-X reference matrix derived from the public NPC Hi-C source data. |
| `data/processed/NPC_chrX/list_muS_unified_260614.RData` | Combined manuscript-scale matrix object for Raw, HiCBZIP-GB/GB(NB), HiCBZIP-N(M), HiCImpute, scHiCluster, Higashi, and Fast-Higashi heatmap and metric summaries. |

## Workflow

Fast path from archived data:

1. Place `list_muS_unified_260614.RData` in `data/processed/NPC_chrX/`.
2. Run `05_summarize_NPC_chrX_manuscript_metrics.R`.
   - Outputs are written under `results/NPC_chrX/manuscript_summaries/` and `results/NPC_chrX/figures/`.

HiCBZIP rerun path:

1. Place `NPC250k_0h_X.mat` in `data/processed/NPC_chrX/`.
2. Run `01_make_NPC_chrX_processed_input.R`.
   - Output: `data/processed/NPC_chrX/data_NPC250k_0h_X_full.RData`.
3. Run `02_run_HiCBZIP_GB_NB_NPC_chrX_all_coverage.R`.
   - Output: `results/NPC_chrX/HiCBZIP_GB_NB/o.bziphic_gbnb.muS.full.AllCoverage.RData`.
4. Run `03_run_HiCBZIP_NM_NPC_chrX_one_coverage.R` for each reported downsampling level.
   - Outputs are written under `results/NPC_chrX/HiCBZIP_NM/`.
5. To rebuild the combined object from method outputs, run `04_build_NPC_chrX_combined_method_object.R`.
   - External benchmark outputs are expected under `results/NPC_chrX/` using the file names shown in that script.
6. Run `05_summarize_NPC_chrX_manuscript_metrics.R`.
   - Outputs are written under `results/NPC_chrX/manuscript_summaries/` and `results/NPC_chrX/figures/`.

The numbered scripts are the primary public entry points for this study. Detailed benchmark-method workflows are kept in `benchmark_methods/` for provenance and are not required when starting from the archived unified object.
