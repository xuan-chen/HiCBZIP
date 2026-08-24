# Run HiCBZIP-GB/GB(NB) for NPC chrX over the manuscript downsampling levels.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "magrittr"))
source_hicbzip_core()

input_file <- Sys.getenv(
  "HICBZIP_REAL1_INPUT_RDATA",
  unset = path_here("data", "processed", "NPC_chrX", "data_NPC250k_0h_X_full.RData")
)
out_dir <- path_here("results", "NPC_chrX", "HiCBZIP_GB_NB")
out_file <- file.path(out_dir, "o.bziphic_gbnb.muS.full.AllCoverage.RData")

require_files(input_file, label = "real-data study 1 processed input")
ensure_dir(out_dir)
load(input_file)
if (!exists("true_muS")) true_muS <- replicate(K, bulk)

coverage_levels <- strsplit(
  Sys.getenv("HICBZIP_REAL1_COVERAGES", unset = "1e-04,2e-04,5e-04,0.001,0.002,0.005,0.01"),
  ",",
  fixed = TRUE
)[[1]]

sim_y_list <- lapply(coverage_levels, function(coverage) {
  set.seed(123456)
  matrix(
    rpois(nrow(true_muS) * ncol(true_muS), lambda = as.numeric(coverage) * true_muS),
    nrow = nrow(true_muS),
    ncol = ncol(true_muS)
  )
})
names(sim_y_list) <- coverage_levels

muS_bzip_gb_nb <- run_BZIP_GB_NB_list(
  sim_y_list,
  as.numeric(coverage_levels),
  r = 1,
  threshold = 3,
  verbose = TRUE
)

save(muS_bzip_gb_nb, coverage_levels, file = out_file)
message("Saved: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
