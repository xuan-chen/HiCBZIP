# Packages
rm(list=ls())
library(tidyverse)
library(cmdstanr)
setwd("~/schic/Output/output_250628_sim_HBA_3X10X0.X_50M_55M_md5M")
source("~/schic/Functions/Functions.R")
load(file = "~/schic/data/Simulation_snm3Cseq_human_brian_astrocytes_res_50k_chr1_50M_55M_md5M_K3X10_coverage_0.1.RData")


set.seed(123456)
args <- commandArgs(trailingOnly = TRUE)
print(args)
coverage <- as.numeric(args[1]) 
sim_y = matrix(rpois(nrow(true_muS)*ncol(true_muS), lambda = coverage*true_muS),nrow(true_muS),ncol(true_muS))


N = nrow(sim_y)
K = ncol(sim_y)
n = length(bin_range)

file <- file.path("~/schic/Functions/BHZIP_nocov_250629.stan")
mod <- cmdstan_model(file)
dt_stan_input = list(N = N,
                     K = K,
                     lambda = rep(coverage, K),
                     # P_mu = P,
                     # P_pi = P,
                     b_i = rep(0,N),
                     # X_mu = X_stacked[,-1],
                     # X_pi = X_stacked[,-1],
                     Y = sim_y)
stan_fit_full = mod$sample(data = dt_stan_input,
                           chains = 2, parallel_chains = 2, 
                           iter_warmup = 250, iter_sampling = 250, save_warmup = T,
                           thin = 1, refresh = 5)
draws_arr <- stan_fit_full$draws() # or format="array"
str(draws_arr)

save(stan_fit_full, draws_arr, file = paste0("CMDSTAN_sim_HBA_3X10X", coverage, "_50M_55M_md5M_nocov_250628.RData"))

muS_indices <- grep("^muS", dimnames(draws_arr)$variable)
muS_names <- dimnames(draws_arr)$variable[muS_indices]
muS_samples <- draws_arr[, , muS_indices]
impute.cmdstan = apply(muS_samples,3, mean) %>% matrix(., ncol = K)

save(impute.cmdstan, file = paste0("Imputation_CMDSTAN_sim_HBA_3X10X", coverage, "_50M_55M_md5M_nocov_250628.RData"))