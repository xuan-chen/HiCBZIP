# Real-data study 1: downsampled NPC chrX

Workflow:

1. Place `NPC250k_0h_X.mat` in `data/processed/NPC_chrX/`.
2. Run `01_make_NPC_chrX_processed_input.R`.
3. Run `02_run_HiCBZIP_GB_NB_NPC_chrX_all_coverage.R`.
4. Run `03_run_HiCBZIP_NM_NPC_chrX_one_coverage.R` for each reported coverage.
5. Place processed external benchmark outputs in the paths documented in `../data/README.md`.
6. Run `04_summarize_NPC_chrX_manuscript_metrics.R`.

The numbered scripts are the primary entry points for reproducing this study. Additional notebooks in this directory preserve the final manuscript summary workflows.
