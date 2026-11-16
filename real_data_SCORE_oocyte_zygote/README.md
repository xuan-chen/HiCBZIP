# Real-Data Study 2: SCORE Oocyte-To-Zygote Benchmark

This folder contains the public workflow for the mouse oocyte-to-zygote downstream benchmark. The analysis uses SCORE-formatted 1 Mb `mm10` input data and compares external-imputation pipelines against integrated embedding methods.

## Fixed-Cohort HiCBZIP Fits

1. Provide the SCORE raw pair files in `oocyte_zygote_mm10/1M/`, the anchor file
   `mm10.genome_split_1M`, and `oocyte_zygote_ref` under
   `data/processed/SCORE_oocyte_zygote/source_inputs/` (or set
   `HICBZIP_SCORE_SOURCE`). The [raw SCORE QC table](../data/processed/SCORE_oocyte_zygote/source_inputs/oocyte_zygote_raw_pairs_filtered_ref)
   is included and used by default, even when the raw inputs are stored
   elsewhere. It records the 150 cells that passed and the two excluded cells.
   Set `HICBZIP_SCORE_RAW_QC` only to override its path.
2. Run `Rscript real_data_SCORE_oocyte_zygote/01_build_SCORE_HiCBZIP_inputs.R`
   from the release root. This freezes the cohort, sums intrachromosomal
   contacts across all chromosomes including the diagonal, and computes
   whole-cell exposure. It then fits GB(NB) and saves observed counts,
   hyperparameters, and unrounded estimates in `gbnb_fits/`.
   Set `HICBZIP_SCORE_PREPARE_ONLY=1` to write just the cohort and exposure
   metadata. Set `HICBZIP_GBNB_CHROMOSOMES` to a comma-separated list to
   limit fitting without changing the cohort or exposure calculation.
3. On Linux with CmdStan configured, run
   `Rscript real_data_SCORE_oocyte_zygote/build_SCORE_HiCBZIP_NM_inputs.R`.
   This maps the saved GB(NB) parameters to delta-method priors and fits N(M)
   using one chain per locus pair (500 warmup and 500 sampling iterations).
   It saves unrounded outputs in `nm_fits/` and can resume from checkpoints.
   Use `HICBZIP_NM_CHROMOSOMES` and `HICBZIP_NM_WORKERS` to limit fitting.

Both scripts default to `results/SCORE_oocyte_zygote/fixed150/`. Set
`HICBZIP_SCORE_FIXED150_OUTPUT` consistently for both stages to change it.

## Manuscript SCORE Benchmark

1. Set `HICBZIP_SCORE_HICIMPUTE_PAIRS` to the preserved HiCImpute pair files
   (one file per fixed-cohort cell, including the four empty outputs).
   Run `Rscript real_data_SCORE_oocyte_zygote/04_prepare_fixed150_depth_matched.R`.
   The script writes continuous InnerProduct inputs and largest-remainder
   integer SnapATAC inputs in `depth_matched_inputs/pair_files/`. It keeps
   150 cells for both HiCBZIP variants and writes a separate reference for
   the 146 nonempty HiCImpute outputs. No raw values substitute for failures.
2. With the benchmark's SCORE environment active, run
   `bash real_data_SCORE_oocyte_zygote/05_run_fixed150_depth_matched_score.sh`.
   The script builds integer `.scool` inputs and runs ten seeds (2026--2035)
   per workflow. Both steps accept `gbnb`, `nm`, or `schicimpute` as
   optional arguments to limit processing. Use a fresh output directory
   when changing the fit or input files, so cached SCORE runs are not reused.
3. Render `06_summarize_fixed150_manuscript.Rmd` in the SCORE analysis
   environment. Set `HICBZIP_SCORE_FIXED150_OUTPUT` to the output root and
   `HICBZIP_SCORE_REPORTED_RESULTS` to the saved Raw and integrated SCORE
   results root (containing `server_json_runs_nm/`,
   `oocyte_zygote_schicluster_10runs/`, `oocyte_zygote_fast_higashi_10runs/`,
   and `oocyte_zygote_higashi_10runs/`). Set `HICBZIP_SCORE_PCA_REFERENCE`
   to the original depth-filtered SCORE reference used for PCA labels when
   reproducing from the archived benchmark. The notebook expects the original
   JSON and embedding layout and writes the 110-run metric summary and PCA
   figures to `manuscript_summaries/`. Use only trusted SCORE embedding files.

The main analysis uses the saved Raw and integrated runs. The scripts here
do not refit HiCImpute, rerun the integrated comparators, or reproduce the
native/direct-rounding sensitivity analysis.

**SnapATAC weighting.** The benchmark used SCORE 0.1.0 and SnapATAC2 2.8.0
(see the [environment record](../environment/conda_list_server_score.txt)).
In the preserved JinLabBioinfo/SCORE checkout, based on commit
`146738d992747ea4d63052d4a8d6fb6623e56c3e`, passing `--snapatac_no_idf`
passes `feature_weights=None` to `snap.tl.spectral()`, which applies IDF
weighting. Omitting the flag supplies all-ones weights. The release runner
retains the flag to reproduce the reported IDF-weighted analysis. Verify
this behavior before substituting another SCORE implementation.

