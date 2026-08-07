# Packages
rm(list=ls())
library(tidyverse)
library(cmdstanr)
setwd("~/schic/Output/output_250628_sim_HBA_3X10X0.X_50M_55M_md5M_GM")
source("~/schic/Functions/Functions.R")
load(file = "~/schic/data/Simulation_snm3Cseq_human_brian_astrocytes_res_50k_chr1_50M_55M_md5M_K3X10_coverage_0.1.RData")


set.seed(123456)
args <- commandArgs(trailingOnly = TRUE)
coverage <- as.numeric(args[1]) 
sim_y = matrix(rpois(nrow(true_muS)*ncol(true_muS), lambda = coverage*true_muS),nrow(true_muS),ncol(true_muS))


N = nrow(sim_y)
K = ncol(sim_y)
n = length(bin_range)
tau_a = 1

gaussian_kernel <- function(size, sigma) {
  center <- (size - 1) / 2
  x <- -center:center
  kernel <- exp(-(x^2) / (2 * sigma^2))
  kernel <- kernel %o% kernel # Outer product for 2D kernel
  kernel <- kernel / sum(kernel) # Normalize
  return(kernel)
}

# Convolution function using base R
convolve2d <- function(matrix, kernel) {
  pad_size <- floor(nrow(kernel) / 2)
  padded_matrix <- matrix(0, nrow(matrix) + 2 * pad_size, ncol(matrix) + 2 * pad_size)
  padded_matrix[(pad_size + 1):(nrow(padded_matrix) - pad_size), 
                (pad_size + 1):(ncol(padded_matrix) - pad_size)] <- matrix
  
  result <- matrix(0, nrow(matrix), ncol(matrix))
  for (i in 1:nrow(result)) {
    for (j in 1:ncol(result)) {
      sub_matrix <- padded_matrix[i:(i + 2 * pad_size), j:(j + 2 * pad_size)]
      result[i, j] <- sum(sub_matrix * kernel)
    }
  }
  return(result)
}
input_matrix <- matrix_long_to_matrix2D_offdiag(rowMeans(sim_y))
kernel <- gaussian_kernel(size = 5, sigma = 1)  # 5x5 kernel, sigma = 1
smoothed_matrix <- convolve2d(input_matrix, kernel)
m_a = log(matrix2D_to_matrix_long(smoothed_matrix, F)+min(smoothed_matrix[smoothed_matrix>0]))



file <- file.path("~/schic/Functions/BHZIP_GM_nocov_250704.stan")
mod <- cmdstan_model(file)
dt_stan_input = list(N = N,
                     K = K,
                     lambda = rep(coverage, K),
                     # P_mu = P,
                     # P_pi = P,
                     b_i = rep(0,N),
                     theta = m_a,
                     tau = tau_a,
                     # X_mu = X_stacked[,-1],
                     # X_pi = X_stacked[,-1],
                     Y = sim_y)
stan_fit_full = mod$sample(data = dt_stan_input, seed = 123456,
                           chains = 2, parallel_chains = 2, 
                           iter_warmup = 250, iter_sampling = 250, save_warmup = T,
                           thin = 1, refresh = 5)
draws_arr <- stan_fit_full$draws() # or format="array"
str(draws_arr)

muS_indices <- grep("^muS", dimnames(draws_arr)$variable)
muS_names <- dimnames(draws_arr)$variable[muS_indices]
muS_samples <- draws_arr[, , muS_indices]
impute.cmdstan = apply(muS_samples,3, mean) %>% matrix(., ncol = K)

save(impute.cmdstan, file = paste0("Imputation_CMDSTAN_sim_HBA_3X10X", coverage, "_50M_55M_md5M_GM_nocov_250704.RData"))