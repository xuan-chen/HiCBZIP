# HiCBZIP clean code draft

This folder is a public-release draft for the HiCBZIP manuscript repository. It is organized around the manuscript workflow:

1. Core HiCBZIP method code.
2. Simulation study.
3. Real-data study 1: downsampled NPC chrX.
4. Real-data study 2: SCORE mouse oocyte-to-zygote benchmark.

The original working code in `../code/` is intentionally left untouched.

## Main entry points

- `WORKFLOW.md`: ordered description of all manuscript steps.
- `HiCBZIP/`: core HiCBZIP implementation and shared analysis helpers.
- `simulations/01_run_HiCBZIP_GB_NB_simulation.R`: simulation HiCBZIP-GB/GB(NB) runner.
- `simulations/02_prepare_simulation_manuscript_summaries.R`: simulation manuscript summary workflow.
- `real_data_study_1/01_make_NPC_chrX_processed_input.R`: NPC chrX processed input builder.
- `real_data_study_1/02_run_HiCBZIP_GB_NB_NPC_chrX_all_coverage.R`: NPC chrX HiCBZIP-GB/GB(NB) runner.
- `real_data_study_1/03_run_HiCBZIP_NM_NPC_chrX_one_coverage.R`: NPC chrX HiCBZIP-N(M) runner.
- `real_data_study_1/04_summarize_NPC_chrX_manuscript_metrics.R`: NPC chrX manuscript summary workflow.
- `real_data_study_2/01_build_SCORE_HiCBZIP_inputs.R`: SCORE Raw/HiCBZIP `.scool` builder.
- `real_data_study_2/02_run_SCORE_embeddings.R`: SCORE embedding runs.
- `real_data_study_2/03_summarize_SCORE_manuscript_metrics.R`: SCORE manuscript metric summary.

## Data policy

Generated figure files are not intended to be committed. Small processed examples can live in `data/`; large processed inputs and result folders should be archived externally and linked from the manuscript and release notes.

See `data/README.md` and `../missing_inputs.md` before copying this folder into the public GitHub repository.
