# Real-data study 2: SCORE benchmark

Ordered public workflow:

1. Place SCORE processed raw inputs under `data/processed/SCORE_oocyte_zygote/`.
2. Run `01_build_SCORE_HiCBZIP_inputs.R`.
3. Add the processed scHiCImpute `.scool` input.
4. Run `02_run_SCORE_embeddings.R`.
5. Add processed integrated-method JSON result folders for scHiCluster, Higashi, and Fast-Higashi.
6. Run `03_summarize_SCORE_manuscript_metrics.R`.

Included source workflows:

- `build_SCORE_HiCBZIP_GB_inputs.Rmd`: builds Raw and HiCBZIP-GB/GB(NB) SCORE-ready pair files and `.scool` inputs.
- `build_SCORE_HiCBZIP_NM_inputs.R`: builds HiCBZIP-N(M) SCORE-ready pair files and `.scool` inputs.
- `run_SCORE_innerproduct_four_inputs_10runs.R`: runs SCORE InnerProduct over Raw, HiCBZIP-GB/GB(NB), HiCBZIP-N(M), and scHiCImpute.
- `run_SCORE_snapatac_four_inputs_10runs.R`: runs SCORE SnapATAC/no-IDF over the same four inputs.
- `summarize_SCORE_final_metrics.Rmd`: final manuscript metric summary.
