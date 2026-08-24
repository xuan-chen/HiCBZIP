# Run HiCBZIP-GB/GB(NB) on the processed simulation inputs used in the manuscript.
#
# Expected input:
#   data/processed/simulation/input/Simulation_snm3Cseq_human_brian_astrocytes_50k_chr*_K3X10.RData
# Each RData file must contain `true_muS`, an N-pair by K-cell ground-truth matrix.
#
# Output:
#   results/simulation/HiCBZIP_GB_NB/BZIP_N_GB_muS_allchr_allcov.rds

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "magrittr"))
source_hicbzip_core()

input_dir <- Sys.getenv(
  "HICBZIP_SIM_INPUT_DIR",
  unset = path_here("data", "processed", "simulation", "input")
)
output_dir <- Sys.getenv(
  "HICBZIP_SIM_OUTPUT_DIR",
  unset = path_here("results", "simulation", "HiCBZIP_GB_NB")
)
ensure_dir(output_dir)

lambda_list <- as.numeric(strsplit(
  Sys.getenv("HICBZIP_SIM_COVERAGES", unset = "0.01,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.7,1"),
  ",",
  fixed = TRUE
)[[1]])

list_sim_files <- function(dir) {
  patt <- "^Simulation_snm3Cseq_human_brian_astrocytes_50k_chr[0-9XY]+_.*_K3X10\\.RData$"
  list.files(dir, pattern = patt, full.names = TRUE)
}

extract_chr <- function(fname) sub(".*chr([0-9XY]+).*", "\\1", basename(fname))

infer_pair_geometry <- function(n_pairs) {
  n_diag <- as.integer(round((sqrt(1 + 8 * n_pairs) - 1) / 2))
  if (n_diag * (n_diag + 1) / 2 == n_pairs) {
    return(list(n_bins = n_diag, include_diag = TRUE))
  }
  n_offdiag <- as.integer(round((1 + sqrt(1 + 8 * n_pairs)) / 2))
  if (n_offdiag * (n_offdiag - 1) / 2 == n_pairs) {
    return(list(n_bins = n_offdiag, include_diag = FALSE))
  }
  stop("Cannot infer number of bins from ", n_pairs, " pair entries.", call. = FALSE)
}

sim_files <- list_sim_files(input_dir)
if (length(sim_files) == 0) {
  stop("No processed simulation RData files found in: ", input_dir, call. = FALSE)
}

muS_all <- list()
for (f in sim_files) {
  env <- new.env(parent = emptyenv())
  load(f, envir = env)
  if (!exists("true_muS", envir = env)) stop("`true_muS` not found in ", f, call. = FALSE)

  true_muS <- as.matrix(env$true_muS)
  chr <- extract_chr(f)
  pair_geometry <- infer_pair_geometry(nrow(true_muS))
  n_bins <- pair_geometry$n_bins
  include_diag <- pair_geometry$include_diag
  list_neighbor_id <- get_neighbor_id(n_bins, r = 1, include_diag = include_diag)

  muS_all[[chr]] <- setNames(vector("list", length(lambda_list)), as.character(lambda_list))
  for (coverage in lambda_list) {
    set.seed(123456)
    sim_y <- matrix(
      rpois(nrow(true_muS) * ncol(true_muS), lambda = coverage * as.numeric(true_muS)),
      nrow = nrow(true_muS),
      ncol = ncol(true_muS)
    )
    message("Running HiCBZIP-GB/GB(NB): chr", chr, ", coverage=", coverage)
    muS_all[[chr]][[as.character(coverage)]] <- run_BZIP_GB_NB(
      sim_y = sim_y,
      coverage = coverage,
      r = 1,
      threshold = 3,
      B = 10,
      include_diag = include_diag
    )
  }
}

out_file <- file.path(output_dir, "BZIP_N_GB_muS_allchr_allcov.rds")
saveRDS(muS_all, out_file)
message("Saved: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
