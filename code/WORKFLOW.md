# HiCBZIP manuscript workflow

This clean draft is organized around the manuscript rather than the chronology of the working scripts.

## Core method

1. `HiCBZIP/BHZIP_EB.R`
   - Implements the empirical-Bayes HiCBZIP-GB and HiCBZIP-GB(NB) estimators.
   - Provides `run_BZIP_GB_NB()` and `run_BZIP_GB_NB_list()`.
2. `HiCBZIP/BHZIP_match_normal.stan`
   - CmdStan model for the matched Normal-prior HiCBZIP-N(M) estimator.
3. `HiCBZIP/analysis_helpers.R`
   - Matrix conversion, evaluation, heatmap, SCC, TopDom, and insulation-score helpers used by the manuscript workflows.

## Simulation study

Manuscript target: biologically realistic simulation across chromosome windows and downsampling levels.

Ordered steps:

1. Put processed simulation RData files in `data/processed/simulation/input/`.
   - Each file should contain `true_muS`.
   - Raw simulation regeneration notebooks are intentionally not part of the public workflow.
2. Run `simulations/01_run_HiCBZIP_GB_NB_simulation.R`.
   - Generates `results/simulation/HiCBZIP_GB_NB/BZIP_N_GB_muS_allchr_allcov.rds`.
3. Add processed outputs for methods that are not regenerated in this repo.
   - Examples: HiCBZIP-N(GS), HiCBZIP-N(M), HiCImpute, scHiCluster, Higashi, Fast-Higashi.
4. Run `simulations/02_prepare_simulation_manuscript_summaries.R`.
   - Renders SMSE, SCC, insulation-score, clustering, and heatmap summary workflows.

Important gap found: the old simulation N(M)/GM runner references `BHZIP_GM_nocov_250704.stan`, but that Stan file is not present in the candidate code folder. If the final simulation uses that exact model, add the Stan file before release or replace the runner with the matched-normal model in `HiCBZIP/BHZIP_match_normal.stan`.

## Real-data study 1: NPC chrX

Manuscript target: downsampled NPC bulk Hi-C chromosome X recovery benchmark.

Ordered steps:

1. Put `NPC250k_0h_X.mat` in `data/processed/real_data_study_1/`.
2. Run `real_data_study_1/01_make_NPC_chrX_processed_input.R`.
   - Generates `data/processed/real_data_study_1/data_NPC250k_0h_X_full.RData`.
3. Run `real_data_study_1/02_run_HiCBZIP_GB_NB_NPC_chrX_all_coverage.R`.
   - Generates GB/NB imputation results for all downsampling levels.
4. Run `real_data_study_1/03_run_HiCBZIP_NM_NPC_chrX_one_coverage.R` for each reported coverage.
   - This is computationally heavier and uses CmdStan.
5. Add processed external benchmark outputs for scHiCluster, HiCImpute, Higashi, and Fast-Higashi.
6. Run `real_data_study_1/04_summarize_NPC_chrX_manuscript_metrics.R`.
   - Renders the SCC, insulation-score, and heatmap workflows.

## Real-data study 2: SCORE mouse oocyte-to-zygote

Manuscript target: downstream cell-state separation metrics from SCORE.

Ordered steps:

1. Put SCORE processed inputs under `data/processed/real_data_study_2/`.
   - `oocyte_zygote_mm10/1M/`
   - `mm10.genome_split_1M`
   - `oocyte_zygote_ref`
2. Run `real_data_study_2/01_build_SCORE_HiCBZIP_inputs.R`.
   - Builds Raw, HiCBZIP-GB/GB(NB), and HiCBZIP-N(M) `.scool` inputs.
3. Add the processed scHiCImpute `.scool` input, because that external method is documented rather than regenerated here.
4. Run `real_data_study_2/02_run_SCORE_embeddings.R`.
   - Runs InnerProduct and SnapATAC/no-IDF over Raw, HiCBZIP-GB/GB(NB), HiCBZIP-N(M), and scHiCImpute.
5. Add processed metric JSON folders for integrated external methods.
   - scHiCluster, Higashi, and Fast-Higashi.
6. Run `real_data_study_2/03_summarize_SCORE_manuscript_metrics.R`.
   - Renders the final SCORE metric summary workflow.

## What should not be included

- Generated figure PNGs.
- Local `.Rhistory`, HTML notebooks, logs, or cache files.
- Full external method working directories unless a reviewer specifically requires them.
