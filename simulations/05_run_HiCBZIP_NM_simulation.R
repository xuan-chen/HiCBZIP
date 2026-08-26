# Run HiCBZIP-N(M) on the processed simulation inputs used in the manuscript.
#
# This is the cleaned version of the original
# `CMDSTAN_mclapply_BHZIP_matchN_250921.R` workflow, generalized from one
# hard-coded chromosome window to all processed simulation windows.
#
# Expected input:
#   data/processed/simulation/input/Simulation_snm3Cseq_human_brain_astrocytes_50k_chr*_K3X10.RData
# Each RData file must contain `true_muS`, an N-pair by K-cell ground-truth matrix.
#
# Output:
#   results/simulation/HiCBZIP_NM/HiCBZIP_NM_chr<chr>_coverage_<coverage>.RData
#
# Optional environment variables:
#   HICBZIP_SIM_INPUT_DIR       Override processed simulation input directory.
#   HICBZIP_SIM_OUTPUT_DIR      Override output directory.
#   HICBZIP_SIM_COVERAGES       Comma-separated coverages; default manuscript grid.
#   HICBZIP_SIM_TARGET_CHR      Run one chromosome only, e.g. "1" or "chr1".
#   HICBZIP_STAN_WARMUP         CmdStan warmup iterations; default 500.
#   HICBZIP_STAN_SAMPLING       CmdStan sampling iterations; default 500.
#   HICBZIP_STAN_CHAINS         CmdStan chains per row; default 1.
#   HICBZIP_NM_WORKERS          Parallel row workers on Linux; default cores - 1.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "magrittr", "cmdstanr", "parallel"))
source_hicbzip_core()
suppressPackageStartupMessages(library(parallel))

input_dir <- Sys.getenv(
  "HICBZIP_SIM_INPUT_DIR",
  unset = path_here("data", "processed", "simulation", "input")
)
output_dir <- Sys.getenv(
  "HICBZIP_SIM_OUTPUT_DIR",
  unset = path_here("results", "simulation", "HiCBZIP_NM")
)
ensure_dir(output_dir)

stan_file <- path_here("HiCBZIP", "BHZIP_match_normal.stan")
require_files(stan_file, label = "HiCBZIP-N(M) Stan model")

lambda_list <- as.numeric(strsplit(
  Sys.getenv("HICBZIP_SIM_COVERAGES", unset = "0.01,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.7,1"),
  ",",
  fixed = TRUE
)[[1]])

iter_warmup <- as.integer(Sys.getenv("HICBZIP_STAN_WARMUP", unset = "500"))
iter_sampling <- as.integer(Sys.getenv("HICBZIP_STAN_SAMPLING", unset = "500"))
chains <- as.integer(Sys.getenv("HICBZIP_STAN_CHAINS", unset = "1"))
parallel_chains <- chains
threads_per_chain <- 1

worker_default <- max(1, parallel::detectCores() - 1)
workers <- as.integer(Sys.getenv("HICBZIP_NM_WORKERS", unset = as.character(worker_default)))

list_sim_files <- function(dir) {
  patt <- "^Simulation_snm3Cseq_human_brain_astrocytes_50k_chr[0-9XY]+_.*_K3X10\\.RData$"
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

fit_one_row_nm <- function(i, sim_y, coverage, lambda_by_cell, neighbor_id, stan_mod,
                           chains, parallel_chains, iter_warmup, iter_sampling) {
  Y_sim <- sim_y[i, ]

  # Original 250921 rule: borrow neighboring rows only when nonzero count <= 2.
  Y_input <- if (sum(Y_sim > 0) > 2) {
    Y_sim
  } else {
    as.vector(sim_y[neighbor_id[[i]], , drop = FALSE])
  }

  prior <- get_EBE_Matched_N_prior(
    Y_input,
    c_plus_d = 10,
    lambda = coverage,
    fix_negative_w = TRUE
  )

  fit <- run_cmdstan_sample(
    stan_mod,
    list(
      data = list(
        N = ncol(sim_y),
        Y = as.integer(Y_sim),
        lambda = lambda_by_cell,
        a_norm = prior$a_norm,
        sigma_mu = prior$sigma_mu,
        b_norm = prior$b_norm,
        sigma_pi = prior$sigma_pi
      ),
      seed = 12345 + i,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      save_warmup = FALSE,
      thin = 1,
      refresh = 0,
      threads_per_chain = threads_per_chain
    ),
    silent_stan = TRUE
  )

  fit$summary(variables = "mu_tilde")$mean
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

Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

mod <- cmdstanr::cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

for (f in sim_files) {
  env <- new.env(parent = emptyenv())
  load(f, envir = env)
  if (!exists("true_muS", envir = env)) stop("`true_muS` not found in ", f, call. = FALSE)

  true_muS <- as.matrix(env$true_muS)
  chr <- extract_chr(f)
  geometry <- infer_pair_geometry(nrow(true_muS))
  neighbor_id <- get_neighbor_id(geometry$n_bins, r = 2, include_diag = geometry$include_diag)

  for (coverage in lambda_list) {
    out_file <- file.path(
      output_dir,
      paste0("HiCBZIP_NM_chr", chr, "_coverage_", format_coverage(coverage), ".RData")
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
    lambda_by_cell <- rep(coverage, ncol(sim_y))

    message("Running HiCBZIP-N(M): chr", chr, ", coverage=", coverage, ", workers=", workers)
    muS_rows <- parallel::mclapply(
      seq_len(nrow(sim_y)),
      fit_one_row_nm,
      sim_y = sim_y,
      coverage = coverage,
      lambda_by_cell = lambda_by_cell,
      neighbor_id = neighbor_id,
      stan_mod = mod,
      chains = chains,
      parallel_chains = parallel_chains,
      iter_warmup = iter_warmup,
      iter_sampling = iter_sampling,
      mc.cores = workers,
      mc.preschedule = TRUE
    )

    res_list <- muS_rows
    impute.cmdstan <- do.call(rbind, muS_rows)

    save(
      res_list,
      impute.cmdstan,
      sim_y,
      coverage,
      chr,
      geometry,
      iter_warmup,
      iter_sampling,
      chains,
      workers,
      file = out_file
    )
    message("Saved: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
  }
}
