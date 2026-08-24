# Build SCORE-ready .scool inputs for the SCORE oocyte-to-zygote study.
#
# This runs the HiCBZIP-GB/GB(NB) and HiCBZIP-N(M) input-building workflows.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))

required_inputs <- c(
  path_here("data", "processed", "SCORE_oocyte_zygote", "oocyte_zygote_mm10", "1M"),
  path_here("data", "processed", "SCORE_oocyte_zygote", "mm10.genome_split_1M"),
  path_here("data", "processed", "SCORE_oocyte_zygote", "oocyte_zygote_ref"),
  path_here("HiCBZIP", "BHZIP_match_normal.stan")
)
require_dirs(required_inputs[1], label = "SCORE raw pair directory")
require_files(required_inputs[-1], label = "SCORE processed input")

render_rmd(path_here("real_data_SCORE_oocyte_zygote", "build_SCORE_HiCBZIP_GB_inputs.Rmd"),
           output_dir = path_here("results", "SCORE_oocyte_zygote", "input_build_reports"))

source(path_here("real_data_SCORE_oocyte_zygote", "build_SCORE_HiCBZIP_NM_inputs.R"))
