# Data layout

Small example data can be committed here. Large processed data should be archived externally and downloaded or copied into the same paths before running the workflow.

Expected processed inputs:

| Study | Path | Contents |
| --- | --- | --- |
| Simulation | `data/processed/simulation/input/` | Processed simulation RData files with `true_muS` |
| Simulation | `data/processed/simulation/muS_combined_plot_260124.RData` | Combined processed method outputs for SMSE/SCC/IS summaries |
| Simulation | `data/processed/simulation/muS_combined_plot_260209.RData` | Combined processed method outputs for clustering summaries |
| Simulation | `data/processed/simulation/higashi_unified_allchr_260208.RData` | Processed Higashi simulation output |
| Simulation | `data/processed/simulation/muS_combined_plot_both_fasthigashi_260613.RData` | Processed heatmap-panel input |
| Real 1 | `data/processed/real_data_study_1/NPC250k_0h_X.mat` | NPC chrX bulk contact matrix |
| Real 1 | `data/processed/real_data_study_1/list_muS_unified_260614.RData` | Combined NPC method output for heatmap panels |
| Real 2 | `data/processed/real_data_study_2/oocyte_zygote_mm10/1M/` | SCORE raw pair files |
| Real 2 | `data/processed/real_data_study_2/mm10.genome_split_1M` | SCORE anchor file |
| Real 2 | `data/processed/real_data_study_2/oocyte_zygote_ref` | SCORE cell metadata/reference |

Do not commit large `.scool`, `.RData`, or raw pair archives unless the target repository policy allows large files or Git LFS.
