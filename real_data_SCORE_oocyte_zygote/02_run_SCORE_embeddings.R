# Run SCORE InnerProduct and SnapATAC/no-IDF for the four external-imputation inputs.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))

required_inputs <- c(
  path_here("data", "bhzip_score_compare_from_pairs", "oocyte_zygote_raw_1M.scool"),
  path_here("data", "bhzip_score_compare_from_pairs", "oocyte_zygote_bhzip_1M.scool"),
  path_here("data", "bhzip_score_compare_from_pairs_nm", "oocyte_zygote_bhzip_nm_1M.scool"),
  path_here("data", "schicimpute_score_compare_from_pairs", "oocyte_zygote_schicimpute_1M.scool")
)
require_files(required_inputs, label = "SCORE .scool input")

source(path_here("real_data_SCORE_oocyte_zygote", "run_SCORE_innerproduct_four_inputs_10runs.R"))
source(path_here("real_data_SCORE_oocyte_zygote", "run_SCORE_snapatac_four_inputs_10runs.R"))
