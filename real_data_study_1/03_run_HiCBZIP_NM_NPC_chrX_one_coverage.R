# Run HiCBZIP-N(M) for NPC chrX at one downsampling level.
#
# Usage:
#   Rscript real_data_study_1/03_run_HiCBZIP_NM_NPC_chrX_one_coverage.R 0.01

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "magrittr", "cmdstanr", "parallel"))
source_hicbzip_core()

args <- commandArgs(trailingOnly = TRUE)
coverage <- suppressWarnings(as.numeric(args[1]))
if (is.na(coverage)) coverage <- as.numeric(Sys.getenv("HICBZIP_REAL1_NM_COVERAGE", unset = "0.01"))

input_file <- Sys.getenv(
  "HICBZIP_REAL1_INPUT_RDATA",
  unset = path_here("data", "processed", "real_data_study_1", "data_NPC250k_0h_X_full.RData")
)
stan_file <- Sys.getenv(
  "HICBZIP_NM_STAN_FILE",
  unset = path_here("HiCBZIP", "BHZIP_match_normal.stan")
)
out_dir <- path_here("results", "real_data_study_1", "HiCBZIP_NM")

require_files(c(input_file, stan_file), label = "real-data study 1 N(M) input")
ensure_dir(out_dir)
load(input_file)
if (!exists("true_muS")) true_muS <- replicate(K, bulk)

set.seed(123456)
sim_y <- matrix(
  rpois(nrow(true_muS) * ncol(true_muS), lambda = coverage * true_muS),
  nrow = nrow(true_muS),
  ncol = ncol(true_muS)
)

N <- nrow(sim_y)
K <- ncol(sim_y)
n_bins <- if (exists("n")) n else as.integer((sqrt(1 + 8 * N) - 1) / 2)
list_neighbor_id <- get_neighbor_id(n_bins, r = 1, include_diag = TRUE)
lambda_by_row <- rep(coverage, N)
mod <- cmdstanr::cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

fit_one_row <- function(i) {
  Y_sim <- sim_y[i, ]
  if (sum(Y_sim > 0) > 2) {
    Y_input <- Y_sim
  } else {
    Y_input <- as.vector(sim_y[list_neighbor_id[[i]], , drop = FALSE])
  }
  ebe <- get_EBE_ZNB_Gamma_Beta(Y_input, c_plus_d = 10, lambda = lambda_by_row[i], fix_negative_w = TRUE)

  data_i <- list(
    N = K,
    Y = Y_sim,
    lambda = lambda_by_row[i],
    a_norm = log(ebe$a / ebe$b),
    sigma_mu = sqrt(1 / ebe$a),
    b_norm = ifelse(ebe$c == 0, -1e4, ifelse(ebe$d == 0, 1e4, log(ebe$c / ebe$d))),
    sigma_pi = ifelse(ebe$c == 0 | ebe$d == 0, 1e-4,
                      sqrt((ebe$c + ebe$d)^2 / (ebe$c * ebe$d * (ebe$c + ebe$d + 1))))
  )

  fit <- mod$sample(
    data = data_i,
    seed = 12345 + i,
    chains = 1,
    parallel_chains = 1,
    iter_warmup = 500,
    iter_sampling = 500,
    refresh = 0,
    threads_per_chain = 1,
    output_dir = tempdir()
  )
  list(mu_tilde = fit$summary(variables = "mu_tilde")$mean)
}

workers <- as.integer(Sys.getenv("HICBZIP_REAL1_NM_WORKERS", unset = "1"))
res_list <- parallel::mclapply(seq_len(N), fit_one_row, mc.cores = workers, mc.preschedule = TRUE)

out_file <- file.path(out_dir, paste0("CMDSTAN_mclapply_BHZIP_matchN_NPC250k_0h_X_full_", coverage, ".RData"))
save(res_list, coverage, file = out_file)
message("Saved: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
