# Assemble real-data study 1 manuscript metrics and heatmap panels.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))

required_inputs <- c(
  path_here("data", "processed", "NPC_chrX", "data_NPC250k_0h_X_full.RData"),
  path_here("results", "NPC_chrX", "HiCBZIP_GB_NB", "o.bziphic_gbnb.muS.full.AllCoverage.RData")
)
require_files(required_inputs, label = "real-data study 1 summary input")

output_dir <- path_here("results", "NPC_chrX", "manuscript_summaries")
ensure_dir(output_dir)

render_rmd(path_here("real_data_NPC_chrX", "summarize_NPC_chrX_all_coverage.Rmd"), output_dir = output_dir)
render_rmd(path_here("real_data_NPC_chrX", "make_NPC_chrX_heatmap_panel.Rmd"), output_dir = output_dir)
