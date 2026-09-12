# HiCBZIP Reproducibility Workflow

This document describes the study-by-study workflow used to reproduce the HiCBZIP manuscript analyses from public source files and processed inputs. Large processed inputs and generated outputs are documented here and in `data/README.md`; they are obtained from the cited public sources or associated data archives before running the full workflows.

## Core method

1. `HiCBZIP/BHZIP_EB.R`
   - Implements the empirical-Bayes HiCBZIP-GB and HiCBZIP-GB(NB) estimators.
   - Provides `run_BZIP_GB_NB()` and `run_BZIP_GB_NB_list()`.
2. `HiCBZIP/BHZIP_match_normal.stan`
   - CmdStan model for the matched Normal-prior HiCBZIP-N(M) estimator.
3. `HiCBZIP/analysis_helpers.R`
   - Matrix conversion, evaluation, heatmap, SCC, TopDom, and insulation-score helpers used by the manuscript workflows.

## Simulation study

Purpose: reproduce the biologically informed simulation across chromosome windows and downsampling levels.

Steps:

1. To regenerate the processed simulation inputs from public raw data, place the three GSE130711 source indexed-contact files in `data/raw/simulation/Raw_Data/`: GSM3749700, GSM3750251, and GSM3751478.
   - Run `simulation/00_prepare_simulation_source_contacts.R`.
   - Run SCL using `simulation/01_run_scl_for_simulation.sh`.
   - Run `simulation/02_generate_processed_simulation_inputs.R`.
2. Alternatively, start from archived processed simulation inputs by placing the RData files in `data/processed/simulation/input/`.
   - Each file contains `true_muS`.
3. Run `simulation/03_run_HiCBZIP_GB_NB_simulation.R`.
   - Generates `results/simulation/HiCBZIP_GB_NB/BZIP_N_GB_muS_allchr_allcov.rds`.
4. Run `simulation/04_run_HiCBZIP_NGS_simulation.R`.
   - Generates one HiCBZIP-N(GS) CmdStan output per chromosome and coverage.
5. Run `simulation/05_run_HiCBZIP_NM_simulation.R`.
   - Generates one HiCBZIP-N(M) CmdStan output per chromosome and coverage.
6. Place processed outputs for benchmark methods in `data/processed/simulation/`.
   - Examples: HiCImpute, scHiCluster, Higashi, Fast-Higashi.
   - Study-specific benchmark-method workflows are provided in `simulation/benchmark_methods/`.
7. Run `simulation/06_prepare_simulation_manuscript_summaries.R`.
   - Renders SMSE, SCC, insulation-score, clustering, and heatmap summary workflows.

## NPC chrX recovery study

Purpose: reproduce the downsampled NPC bulk Hi-C chromosome X recovery benchmark.

Steps:

1. For manuscript figure and metric reproduction, put the archived unified NPC object in `data/processed/NPC_chrX/list_muS_unified_260614.RData`.
   - This object contains Raw, HiCBZIP-GB/GB(NB), HiCBZIP-N(M), HiCImpute, scHiCluster, Higashi, and Fast-Higashi matrices across the reported downsampling levels.
   - Run `real_data_NPC_chrX/05_summarize_NPC_chrX_manuscript_metrics.R`.
2. To rerun the HiCBZIP portions, put the processed NPC chrX reference matrix in `data/processed/NPC_chrX/NPC250k_0h_X.mat`.
   - This matrix is derived from the public NPC Hi-C source data from GSE72697, sample GSM1868576.
3. Run `real_data_NPC_chrX/01_make_NPC_chrX_processed_input.R`.
   - Generates `data/processed/NPC_chrX/data_NPC250k_0h_X_full.RData`.
4. Run `real_data_NPC_chrX/02_run_HiCBZIP_GB_NB_NPC_chrX_all_coverage.R`.
   - Generates GB/NB imputation results for all downsampling levels.
5. Run `real_data_NPC_chrX/03_run_HiCBZIP_NM_NPC_chrX_one_coverage.R` for each reported coverage.
   - This is computationally heavier and uses CmdStan.
6. To rebuild the unified NPC object from all method outputs, run `real_data_NPC_chrX/04_build_NPC_chrX_combined_method_object.R`.
   - Detailed benchmark-method workflows are kept in `real_data_NPC_chrX/benchmark_methods/`.
   - Benchmark-method provenance is documented in `benchmarks/README.md`.

## SCORE mouse oocyte-to-zygote study

Purpose: reproduce the downstream cell-state separation benchmark using SCORE.

Steps:

1. Put SCORE processed inputs under `data/processed/SCORE_oocyte_zygote/`.
   - `oocyte_zygote_mm10/1M/`
   - `mm10.genome_split_1M`
   - `oocyte_zygote_ref`
2. Run `real_data_SCORE_oocyte_zygote/01_build_SCORE_HiCBZIP_inputs.R`.
   - Builds Raw, HiCBZIP-GB/GB(NB), and HiCBZIP-N(M) `.scool` inputs.
3. Build or provide the processed scHiCImpute `.scool` input.
   - The study-specific script is `real_data_SCORE_oocyte_zygote/benchmark_methods/01_build_scHiCImpute_SCORE_inputs.R`.
4. Run `real_data_SCORE_oocyte_zygote/02_run_SCORE_embeddings.R`.
   - Runs InnerProduct and SnapATAC/no-IDF over Raw, HiCBZIP-GB/GB(NB), HiCBZIP-N(M), and scHiCImpute.
5. Run or provide processed metric JSON folders for integrated benchmark methods.
   - Study-specific scripts for scHiCluster, Higashi, and Fast-Higashi are in `real_data_SCORE_oocyte_zygote/benchmark_methods/`.
6. Run `real_data_SCORE_oocyte_zygote/03_summarize_SCORE_manuscript_metrics.R`.
   - Renders the final SCORE metric summary workflow.

## Files Excluded From Git

- Generated figure PNGs.
- `.Rhistory`, HTML notebooks, logs, and cache files.
- Full benchmark method working directories.
- Large `.RData`, `.mat`, `.scool`, and raw contact files, unless intentionally tracked with Git LFS.
