# SCORE Oocyte-to-Zygote Benchmark Methods

This folder contains study-specific benchmark-method workflows for real-data study 2. The main real-data-2 workflow in `real_data_SCORE_oocyte_zygote/` focuses on HiCBZIP input generation, SCORE runs, and final manuscript summaries. These scripts document the additional benchmark methods used for comparison.

## Inputs

Run `real_data_SCORE_oocyte_zygote/01_build_SCORE_HiCBZIP_inputs.R` first, or provide equivalent SCORE-ready files at:

```text
data/bhzip_score_compare_from_pairs/oocyte_zygote_raw_1M.scool
data/bhzip_score_compare_from_pairs/oocyte_zygote_bhzip_1M.scool
data/bhzip_score_compare_from_pairs/oocyte_zygote_ref_min_depth_5000.tsv
data/bhzip_score_compare_from_pairs_nm/oocyte_zygote_raw_1M.scool
data/bhzip_score_compare_from_pairs_nm/oocyte_zygote_bhzip_nm_1M.scool
data/bhzip_score_compare_from_pairs_nm/oocyte_zygote_ref_min_depth_5000.tsv
```

The scHiCImpute script also needs the SCORE-formatted raw pair files:

```text
data/processed/SCORE_oocyte_zygote/oocyte_zygote_mm10/1M/
data/processed/SCORE_oocyte_zygote/mm10.genome_split_1M
data/processed/SCORE_oocyte_zygote/oocyte_zygote_ref
```

## Scripts

| Script | Purpose |
| --- | --- |
| `01_build_scHiCImpute_SCORE_inputs.R` | Runs scHiCImpute on the SCORE raw pair files and builds `oocyte_zygote_schicimpute_1M.scool`. |
| `02_run_scHiCluster_SCORE_10runs.R` | Runs SCORE's scHiCluster embedding workflow for 10 seeds. |
| `03_run_Higashi_SCORE_10runs.R` | Runs SCORE's Higashi embedding workflow for 10 seeds. |
| `04_run_FastHigashi_SCORE_10runs.R` | Runs SCORE's Fast-Higashi embedding workflow for 10 seeds. |
| `score_benchmark_helpers.R` | Shared SCORE command, metric-extraction, and summary helpers. |

## Outputs Consumed By The Manuscript Summary

The final real-data-2 summary script reads per-run SCORE JSON outputs from:

```text
results/SCORE_oocyte_zygote/oocyte_zygote_schicimpute_10runs_innerproduct/
results/SCORE_oocyte_zygote/oocyte_zygote_schicimpute_10runs_snapatac_noidf/
results/SCORE_oocyte_zygote/oocyte_zygote_schicluster_10runs/
results/SCORE_oocyte_zygote/oocyte_zygote_fast_higashi_10runs/
results/SCORE_oocyte_zygote/oocyte_zygote_higashi_10runs/
```

The main HiCBZIP workflow writes Raw, HiCBZIP-GB/GB(NB), and HiCBZIP-N(M) SCORE results to the corresponding `results/SCORE_oocyte_zygote/` folders documented in `../README.md`.

Software versions and environment records are documented under `environment/` and `benchmarks/README.md`.
