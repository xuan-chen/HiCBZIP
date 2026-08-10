# Assemble manuscript simulation summaries from processed method outputs.
#
# This script renders the metric, clustering, and heatmap workflows after checking
# the key input files expected by the manuscript notebooks.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))

summary_input <- path_here("data", "processed", "simulation", "muS_combined_plot_260124.RData")
higashi_input <- path_here("data", "processed", "simulation", "higashi_unified_allchr_260208.RData")
cluster_input <- path_here("data", "processed", "simulation", "muS_combined_plot_260209.RData")
output_dir <- path_here("results", "simulation", "manuscript_summaries")

require_files(c(summary_input, higashi_input, cluster_input), label = "simulation processed summary input")
ensure_dir(output_dir)

render_rmd(path_here("simulations", "summarize_simulation_metrics.Rmd"), output_dir = output_dir)
render_rmd(path_here("simulations", "summarize_simulation_clustering.Rmd"), output_dir = output_dir)
render_rmd(path_here("simulations", "make_simulation_heatmap_panels.Rmd"), output_dir = output_dir)
