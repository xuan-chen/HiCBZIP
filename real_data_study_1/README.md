# Real-data study 1: downsampled NPC chrX

Ordered public workflow:

1. Place `NPC250k_0h_X.mat` in `data/processed/real_data_study_1/`.
2. Run `01_make_NPC_chrX_processed_input.R`.
3. Run `02_run_HiCBZIP_GB_NB_NPC_chrX_all_coverage.R`.
4. Run `03_run_HiCBZIP_NM_NPC_chrX_one_coverage.R` for each reported coverage.
5. Add processed external benchmark outputs.
6. Run `04_summarize_NPC_chrX_manuscript_metrics.R`.

The folder keeps the final summary notebooks, but the public workflow should use the numbered scripts above as the reviewer-facing entry points.
