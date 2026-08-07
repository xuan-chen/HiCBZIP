# Correction: Neighboor_threshold > 2 (previously 4, >=5)

# Packages
rm(list = ls())
library(tidyverse)
library(cmdstanr)
library(parallel)

# Workdir & sources
setwd("~/schic/Output/output_251018_sim_HBA_3X10X0.X_chr10_region5M_md5M_NM")
suppressMessages(source("~/schic/Functions/Functions.R"))
suppressMessages(source("~/schic/Functions/BHZIP_EB_250911.R"))

# Data
suppressMessages(load("~/schic/data/Simulation_snm3Cseq_human_brian_astrocytes_50k_chr10_25M_30M_md5M_K3X10.RData"))

# Args
set.seed(123456)
args <- commandArgs(trailingOnly = TRUE)
coverage <- suppressWarnings(as.numeric(args[1]))
if (is.na(coverage)) coverage <- 0.5

# Simulation
sim_y <- matrix(
  rpois(nrow(true_muS) * ncol(true_muS), lambda = coverage * true_muS),
  nrow(true_muS), ncol(true_muS)
)

N <- nrow(sim_y)
K <- ncol(sim_y)
n <- length(bin_range)
n_rep <- 10

# Neighborhood
r <- 1
threshold_neighbor <- 2
list_neighbor_id <- get_neighbor_id(n, r, include_diag = FALSE)

# Compile model
stan_file <- "~/schic/Functions/BHZIP_match_normal_250906.stan"
mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE), quiet = TRUE)

# Hyper-parameters
lambda <- coverage
lambda_by_row <- if (length(lambda) == 1) rep(lambda, N) else as.numeric(lambda)
c_plus_d <- 10
fix_negative_w <- TRUE

# Prevent over-threading
Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

# Parallel setup
total_cores <- parallel::detectCores()
chains_per_dataset <- 1
threads_per_chain  <- 1
workers <- max(1, floor(total_cores / (chains_per_dataset * threads_per_chain)) - 1)
cat("Cores:", total_cores, "| workers:", workers, "| chains/dataset:", chains_per_dataset, "\n")

fit_one_row <- function(i) {
  Y_sim <- sim_y[i, ]
  if (sum(Y_sim > 0) > threshold_neighbor) {
    Y_input <- Y_sim
  } else {
    Y_input <- as.vector(sim_y[list_neighbor_id[[i]], ])
  }
  
  ebe <- get_EBE_ZNB_Gamma_Beta(Y_input, c_plus_d, lambda_by_row[i], fix_negative_w)
  
  a_hat <- ebe$a; b_hat <- ebe$b; c_hat <- ebe$c; d_hat <- ebe$d
  
  a_norm   <- log(a_hat / b_hat)
  sigma_mu <- sqrt(1 / a_hat)
  b_norm   <- ifelse(c_hat == 0, -1e4,
                     ifelse(d_hat == 0,  1e4, log(c_hat / d_hat)))
  sigma_pi <- ifelse(c_hat == 0 | d_hat == 0, 1e-4,
                     sqrt((c_hat + d_hat)^2 / (c_hat * d_hat * (c_hat + d_hat + 1))))
  
  data_i <- list(
    N        = K,
    Y        = Y_sim,
    lambda   = lambda_by_row[i],
    a_norm   = a_norm,
    sigma_mu = sigma_mu,
    b_norm   = b_norm,
    sigma_pi = sigma_pi
  )
  
  fit <- suppressMessages(
    mod$sample(
      data = data_i,
      seed = 12345 + i,
      chains = chains_per_dataset,
      parallel_chains = chains_per_dataset,
      iter_warmup = 500,
      iter_sampling = 500,
      refresh = 0,                     # <-- no progress output
      threads_per_chain = threads_per_chain,
      show_messages = FALSE,
      output_dir = tempdir()           # discard CmdStan CSVs
    )
  )
  
  list(
    priors = c(a_norm = a_norm, sigma_mu = sigma_mu,
               b_norm = b_norm, sigma_pi = sigma_pi),
    mu        = fit$summary(variables = "mu")$mean,
    pi        = fit$summary(variables = "pi")$mean,
    S         = fit$summary(variables = "S")$mean,
    mu_tilde  = fit$summary(variables = "mu_tilde")$mean
  )
}

# Run silently
res_list <- suppressMessages(
  parallel::mclapply(seq_len(N), fit_one_row,
                     mc.cores = workers, mc.preschedule = TRUE)
)

# Save only results
save(res_list,
     file = paste0("Imputation_CMDSTAN_sim_HBA_3X10X",
                   coverage, "_NM_nocov_251018.RData"))
