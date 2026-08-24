# Run HiCBZIP-N(GS) on the processed simulation inputs used in the manuscript.
#
# This is the cleaned version of the original
# `cmdstan_50M_55M_X_BHZIP_GM_nocov_250704.R` workflow, generalized from one
# hard-coded chromosome window to all processed simulation windows.
#
# Expected input:
#   data/processed/simulation/input/Simulation_snm3Cseq_human_brian_astrocytes_50k_chr*_K3X10.RData
# Each RData file must contain `true_muS`, an N-pair by K-cell ground-truth matrix.
#
# Output:
#   results/simulation/HiCBZIP_NGS/HiCBZIP_NGS_chr<chr>_coverage_<coverage>.RData
#
# Optional environment variables:
#   HICBZIP_SIM_INPUT_DIR       Override processed simulation input directory.
#   HICBZIP_SIM_OUTPUT_DIR      Override output directory.
#   HICBZIP_SIM_COVERAGES       Comma-separated coverages; default manuscript grid.
#   HICBZIP_SIM_TARGET_CHR      Run one chromosome only, e.g. "1" or "chr1".
#   HICBZIP_STAN_WARMUP         CmdStan warmup iterations; default 250.
#   HICBZIP_STAN_SAMPLING       CmdStan sampling iterations; default 250.
#   HICBZIP_STAN_CHAINS         CmdStan chains; default 2.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "magrittr", "cmdstanr"))
source_hicbzip_core()

input_dir <- Sys.getenv(
  "HICBZIP_SIM_INPUT_DIR",
  unset = path_here("data", "processed", "simulation", "input")
)
output_dir <- Sys.getenv(
  "HICBZIP_SIM_OUTPUT_DIR",
  unset = path_here("results", "simulation", "HiCBZIP_NGS")
)
ensure_dir(output_dir)

stan_file <- path_here("HiCBZIP", "BHZIP_GM_nocov.stan")
require_files(stan_file, label = "HiCBZIP-N(GS) Stan model")

lambda_list <- as.numeric(strsplit(
  Sys.getenv("HICBZIP_SIM_COVERAGES", unset = "0.01,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.7,1"),
  ",",
  fixed = TRUE
)[[1]])

iter_warmup <- as.integer(Sys.getenv("HICBZIP_STAN_WARMUP", unset = "250"))
iter_sampling <- as.integer(Sys.getenv("HICBZIP_STAN_SAMPLING", unset = "250"))
chains <- as.integer(Sys.getenv("HICBZIP_STAN_CHAINS", unset = "2"))
parallel_chains <- chains

list_sim_files <- function(dir) {
  patt <- "^Simulation_snm3Cseq_human_brian_astrocytes_50k_chr[0-9XY]+_.*_K3X10\\.RData$"
  list.files(dir, pattern = patt, full.names = TRUE)
}

extract_chr <- function(fname) sub(".*chr([0-9XY]+).*", "\\1", basename(fname))

format_coverage <- function(x) gsub("\\.", "p", as.character(x))

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

target_chr <- Sys.getenv("HICBZIP_SIM_TARGET_CHR", unset = "")
if (nzchar(target_chr)) {
  target_chr <- sub("^chr", "", target_chr, ignore.case = TRUE)
  sim_files <- sim_files[vapply(sim_files, extract_chr, character(1)) == target_chr]
  if (length(sim_files) == 0) stop("No simulation file found for target chromosome: ", target_chr, call. = FALSE)
}

for (f in sim_files) {
  env <- new.env(parent = emptyenv())
  load(f, envir = env)
  if (!exists("true_muS", envir = env)) stop("`true_muS` not found in ", f, call. = FALSE)

  true_muS <- as.matrix(env$true_muS)
  chr <- extract_chr(f)
  geometry <- infer_pair_geometry(nrow(true_muS))

  for (coverage in lambda_list) {
    out_file <- file.path(
      output_dir,
      paste0("HiCBZIP_NGS_chr", chr, "_coverage_", format_coverage(coverage), ".RData")
    )
    if (file.exists(out_file)) {
      message("Skipping existing file: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
      next
    }

    set.seed(123456)
    sim_y <- matrix(
      rpois(nrow(true_muS) * ncol(true_muS), lambda = coverage * as.numeric(true_muS)),
      nrow = nrow(true_muS),
      ncol = ncol(true_muS)
    )

    message("Running HiCBZIP-N(GS): chr", chr, ", coverage=", coverage)
    impute.cmdstan <- run_HiCBZIP_NGS(
      sim_y = sim_y,
      coverage = coverage,
      stan_file = stan_file,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      include_diag = geometry$include_diag
    )

    save(
      impute.cmdstan,
      sim_y,
      coverage,
      chr,
      geometry,
      iter_warmup,
      iter_sampling,
      chains,
      file = out_file
    )
    message("Saved: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
  }
}
