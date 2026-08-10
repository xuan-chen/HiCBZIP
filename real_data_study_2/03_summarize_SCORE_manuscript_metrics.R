# Summarize SCORE manuscript metrics for real-data study 2.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))

required_dirs <- c(
  path_here("results", "real_data_study_2", "server_json_runs"),
  path_here("results", "real_data_study_2", "server_json_runs_nm"),
  path_here("results", "real_data_study_2", "oocyte_zygote_schicimpute_10runs_innerproduct"),
  path_here("results", "real_data_study_2", "oocyte_zygote_schicimpute_10runs_snapatac_noidf"),
  path_here("results", "real_data_study_2", "oocyte_zygote_schicluster_10runs"),
  path_here("results", "real_data_study_2", "oocyte_zygote_fast_higashi_10runs"),
  path_here("results", "real_data_study_2", "oocyte_zygote_higashi_10runs")
)
require_dirs(required_dirs, label = "SCORE metric result directory")

render_rmd(path_here("real_data_study_2", "summarize_SCORE_final_metrics.Rmd"),
           output_dir = path_here("results", "real_data_study_2", "manuscript_summaries"))
