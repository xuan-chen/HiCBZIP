# Packages
getwd()
rm(list=ls())
library(tidyverse)
source("Functions/Functions.R")
source("STAN_BHZIP_lambda/BHZIP_EB_250911.R")

data_dir <- "C:/Users/67402/OneDrive - University of Miami/courses/research/HiC/Data/Methylation_Hi-C_data/sim_HBA3_random_chr"

out_dir  <- "Output/output_251026_sim_HBA_ten_chr_combine"

out_rds  <- file.path(out_dir, "BZIP_N_GB_muS_allchr_allcov.rds")

# ===============================================
# BZIP-N(GB): muS_posterior_2 with neighborhood borrowing
# Discover sim files like in HiCImpute
# ===============================================

list_sim_files <- function(dir) {
  patt <- "^Simulation_snm3Cseq_human_brian_astrocytes_50k_chr[0-9XY]+_.*_K3X10\\.RData$"
  list.files(dir, pattern = patt, full.names = TRUE)
}
extract_chr <- function(fname) sub(".*(chr[0-9XY]+).*", "\\1", fname)
infer_n_from_long <- function(N_pairs) as.integer(round((1 + sqrt(1 + 8 * N_pairs)) / 2))

compute_muS_for_coverage <- function(sim_y, coverage, list_neighbor_id, thresh_nonzero = 3) {
  N <- nrow(sim_y); K <- ncol(sim_y)
  muS <- matrix(0, nrow = N, ncol = K)
  lambda_t <- coverage
  for (i in seq_len(N)) {
    Y_sim <- sim_y[i, ]
    if (sum(Y_sim > 0) >= thresh_nonzero) Y_input <- Y_sim
    else                                  Y_input <- as.vector(sim_y[list_neighbor_id[[i]], , drop = FALSE])
    ebe <- get_EBE_ZNB_Gamma_Beta(Y_input, c_plus_d = 10, lambda = lambda_t)
    if (!is.finite(ebe$a) || ebe$a < 0) {
      ebe_all <- get_EBE_ZNB_Gamma_Beta(as.vector(sim_y), c_plus_d = 10, lambda = lambda_t)
      ebe$a <- ebe_all$a; ebe$b <- ebe_all$b
    }
    E_S     <- get_E_S(ebe$a, ebe$b, ebe$c, ebe$d, lambda_t)
    S_post  <- ifelse(Y_sim > 0, 0, E_S)
    mu_post <- (ebe$a + Y_sim) / (ebe$b + lambda_t)
    muS[i, ] <- ifelse(Y_sim == 0, mu_post * (1 - S_post), mu_post)
  }
  muS
}

lambda_list <- c(0.01,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.7,1)
sim_files <- list_sim_files(data_dir)
if (!length(sim_files)) stop("No simulation files found in: ", data_dir)
message("Found ", length(sim_files), " files.")

total_tasks <- length(sim_files) * length(lambda_list)
pb <- txtProgressBar(min = 0, max = total_tasks, style = 3)
done <- 0

muS_all <- list()

for (f in sim_files) {
  chr <- extract_chr(basename(f))
  message("\n=== ", basename(f), " (", chr, ") ===")
  
  e <- new.env(parent = emptyenv()); load(f, envir = e)
  if (!exists("true_muS", envir = e)) { warning("true_muS missing in ", f); next }
  true_muS <- e$true_muS
  N <- nrow(true_muS); K <- ncol(true_muS)
  
  # neighbors per chromosome
  n_bins <- infer_n_from_long(N)
  list_neighbor_id <- get_neighbor_id(n_bins, r = 1, include_diag = FALSE)
  
  muS_all[[chr]] <- setNames(vector("list", length(lambda_list)), as.character(lambda_list))
  
  for (covg in lambda_list) {
    set.seed(123456)
    sim_y <- matrix(rpois(N * K, lambda = covg * as.numeric(true_muS)), nrow = N, ncol = K)
    t0 <- Sys.time()
    muS_hat <- compute_muS_for_coverage(sim_y, covg, list_neighbor_id, thresh_nonzero = 3)
    cat(sprintf("  coverage = %-4s | %.1f sec\n",
                as.character(covg),
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    flush.console()
    muS_all[[chr]][[as.character(covg)]] <- muS_hat
    done <- done + 1; setTxtProgressBar(pb, done)
  }
}

close(pb)
saveRDS(muS_all, out_rds)
cat("\nSaved: ", normalizePath(out_rds, winslash = "/"), "\n")
