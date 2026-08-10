# Core HiCBZIP empirical Bayes and matched-Normal imputation functions.
#
# Dependencies:
# - stats, utils, parallel: base/recommended R packages
# - dplyr and magrittr: used by neighborhood helper code
# - cmdstanr: required only for get_EBE_Matched_N() / get_EBE_Matched_N_cmd()
#
# Example setup:
# library(dplyr)
# library(magrittr)
# library(cmdstanr)

# log
## 26/01/11
#### 1. Fixed max->pmax in ebe$b in get_EBE_ZNB_Gamma_Beta
#### 2. Swith b_vec to b_mat in run_BZIP_GB_NB



# Objective and EBE Estimation Functions (Method of Moments on nonzero data)
objective <- function(params, ybar, svar) {
  a <- params[1]
  rho <- params[2]
  # Avoid problematic values:
  if(a <= 0 || rho <= 0 || rho >= 1) return(1e10)
  
  # Conditional mean and variance of the NB component (given nonzero)
  model_mean <- a * (1 - rho) / (rho * (1 - rho^a))
  model_var <- (a * (1 - rho) + a^2 * (1 - rho)^2) / (rho^2 * (1 - rho^a)) -
    (a^2 * (1 - rho)^2) / (rho^2 * (1 - rho^a)^2)
  
  err1 <- (model_mean - ybar)^2
  err2 <- (model_var - svar)^2
  return(err1 + err2)
}

get_EBE_ZNB_Gamma_Beta = function(Y, c_plus_d = 10, lambda = 1,fix_negative_w = T) {
  prop_zero <- mean(Y == 0)
  y_nonzero <- Y[Y > 0]
  
  ybar_plus <- mean(y_nonzero)
  svar_plus <- var(y_nonzero)
  
  # Issue 1: deal with zero or one nonzero observation
  if(sum(Y>0) == 1) {
    return(list(a = 1e6, b = 1e6*lambda/mean(Y), c = c_plus_d, d = 0, rho = 0.5, w = 1, ybar_plus = ybar_plus, svar_plus = svar_plus))
  }
  if(sum(Y>0) == 0) {
    return(list(a = 1e6, b = 1e6*lambda/1e-8, c = c_plus_d, d = 0, rho = 0.5, w = 1, ybar_plus = ybar_plus, svar_plus = svar_plus))
  }
  
  
  opt_result <- optim(par = c(1.5, 0.5), 
                      fn = objective, 
                      ybar = ybar_plus, 
                      svar = svar_plus,
                      method = "L-BFGS-B",
                      lower = c(0.001, 0.01), upper = c(1000, 0.99))
  a_hat <- opt_result$par[1]
  rho_hat <- opt_result$par[2]
  
  # Recover b using the relationship: rho = b/(b + lambda)  =>  b = lambda * rho/(1 - rho)
  b_hat <- lambda * rho_hat / (1 - rho_hat)
  
  # Estimate w from the zero probability equation:
  # P(Y = 0) = w + (1 - w) * rho^a   where w = c/(c + d)
  w_hat <- (prop_zero - rho_hat^(a_hat)) / (1 - rho_hat^(a_hat))
  c_hat <- c_plus_d * w_hat
  d_hat <- c_plus_d * (1 - w_hat)
  
  # Issue 2: Deal with negative w_hat #
  if (w_hat<0 & fix_negative_w){
    c_hat <- c_plus_d * 0
    d_hat <- c_plus_d * (1 - 0)
    a_hat = max(0,mean(Y)^2/(var(Y)-mean(Y)))
    b_hat = pmax(0,lambda*mean(Y)/(var(Y)-mean(Y)))
    rho_hat = b_hat/(b_hat+lambda)
  }
  
  if (var(Y) <= mean(Y)){
    a_hat = 1e6
    b_hat = a_hat*lambda/mean(Y)
  }
  return(list(a = a_hat, b = b_hat, c = c_hat, d = d_hat, rho = rho_hat, w = w_hat, ybar_plus = ybar_plus, svar_plus = svar_plus))
}


get_E_S = function(a,b,c,d,lambda){
  c / (c + d * (b/(b+lambda))^ a)
}


get_EBE_Matched_N = function(Y, c_plus_d = 10, lambda = 1, fix_negative_w = T,
                             stan_file = "BHZIP_match_normal.stan") {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop("cmdstanr is required for matched-Normal fitting. Install cmdstanr or pass a compiled model to get_EBE_Matched_N_cmd().")
  }
  stan_mod <- cmdstanr::cmdstan_model(stan_file)
  get_EBE_Matched_N_cmd(Y, c_plus_d, lambda, fix_negative_w, stan_mod)
}

get_MSE = function(truth, pred){
  MSE = mean((pred-truth)^2)
  return(MSE)
}

get_SMRE = function(truth, pred){
  MSE = mean((pred-truth)^2)
  var_truth = var(truth)
  return(MSE/var_truth)
}


get_neighbor_id_slow = function(n,r,include_diag = F){
  if (include_diag){
    idx_info = expand.grid(j=1:n,i=1:n) %>%
      filter(j>=i) %>%
      mutate(id = row_number())
  }else{
    idx_info = expand.grid(j=1:n,i=1:n) %>%
      filter(j>i) %>%
      mutate(id = row_number())
  }
  
  lapply(1:nrow(idx_info), function(t){
    i0 = idx_info$i[t]
    j0 = idx_info$j[t]
    idx_select = idx_info %>%
      filter(abs(i0-i)+abs(j0-j)<=r)
    return(idx_select$id)
  })
}

get_neighbor_id <- function(n, r, include_diag = TRUE) {
  
  # Step 1: Create the coordinate-to-ID mapping.
  # This part is kept identical to your original function to ensure the `id` values
  # are exactly the same. This setup is fast and not the bottleneck.
  if (include_diag) {
    idx_info <- expand.grid(j = 1:n, i = 1:n) %>%
      filter(j >= i) %>%
      mutate(id = row_number())
  } else {
    idx_info <- expand.grid(j = 1:n, i = 1:n) %>%
      filter(j > i) %>%
      mutate(id = row_number())
  }
  
  # Step 2: Create a fast lookup matrix (id_map).
  # This matrix allows us to find the id for any (i, j) coordinate in O(1) time.
  # We initialize an n x n matrix with NAs, then fill in the valid IDs.
  id_map <- matrix(NA_integer_, nrow = n, ncol = n)
  id_map[cbind(idx_info$i, idx_info$j)] <- idx_info$id
  
  # Step 3: Pre-calculate the relative offsets for the neighborhood.
  # Instead of checking distances in a loop, we pre-define the "diamond" shape
  # of neighbors based on the Manhattan distance (abs(di) + abs(dj) <= r).
  rel_coords <- expand.grid(di = -r:r, dj = -r:r)
  rel_coords <- rel_coords[abs(rel_coords$di) + abs(rel_coords$dj) <= r, ]
  
  # Step 4: Efficiently find neighbors for each point.
  # We iterate through each point, calculate neighbor coordinates using vectorized
  # addition, and use the id_map for a quick lookup.
  
  coords_matrix <- as.matrix(idx_info[, c("i", "j")])
  num_points <- nrow(coords_matrix)
  result_list <- vector("list", num_points)
  
  for (k in 1:num_points) {
    i0 <- coords_matrix[k, 1]
    j0 <- coords_matrix[k, 2]
    
    # Calculate absolute coordinates of all potential neighbors at once
    neighbor_i <- i0 + rel_coords$di
    neighbor_j <- j0 + rel_coords$dj
    
    # Quickly filter out coordinates that are outside the n x n grid
    valid_mask <- neighbor_i >= 1 & neighbor_i <= n & 
      neighbor_j >= 1 & neighbor_j <= n
    
    # Check if there are any valid neighbors before proceeding
    if (any(valid_mask)) {
      # Use matrix indexing on the id_map to get all neighbor IDs in one step
      neighbor_ids <- id_map[cbind(neighbor_i[valid_mask], neighbor_j[valid_mask])]
      
      # Store the final list of valid (non-NA) IDs
      result_list[[k]] <- neighbor_ids[!is.na(neighbor_ids)]
    } else {
      # If no valid neighbors, store an empty integer vector
      result_list[[k]] <- integer(0)
    }
  }
  
  return(result_list)
}

get_EBE_Matched_N_cmd = function(Y, c_plus_d = 10, lambda = 1, fix_negative_w = T, stan_mod){
  
  ebe = get_EBE_ZNB_Gamma_Beta(Y, c_plus_d, lambda, fix_negative_w)
  a_hat = ebe$a
  b_hat = ebe$b
  c_hat = ebe$c
  d_hat = ebe$d
  
  a_norm <- log(a_hat / b_hat)
  sigma_mu <- sqrt(1 / a_hat)
  b_norm <- ifelse(c_hat ==0, -1e4, ifelse(d_hat ==0, 1e4, log(c_hat / d_hat))) ## make two options
  sigma_pi <- ifelse(c_hat ==0 | d_hat == 0, 1e-4, sqrt((c_hat + d_hat)^2 / (c_hat * d_hat * (c_hat + d_hat + 1))))
  
  
  data_list_normal <- list(
    N = length(Y),
    Y = Y,
    lambda = lambda,
    a_norm = a_norm,
    sigma_mu = sigma_mu,
    b_norm = b_norm,
    sigma_pi = sigma_pi
  )
  
  stan_fit_full = stan_mod$sample(data = data_list_normal,seed = 12345,
                             chains = 2, parallel_chains = 2, 
                             iter_warmup = 500, iter_sampling = 500, save_warmup = T,
                             thin = 1, refresh = 0)
  
  muS_posterior <- stan_fit_full$summary(variables = "mu_tilde")$mean
  mu_posterior  <- stan_fit_full$summary(variables = "mu")$mean
  pi_posterior  <- stan_fit_full$summary(variables = "pi")$mean
  S_posterior   <- stan_fit_full$summary(variables = "S")$mean
  
  ebe_normal = 
    list(
      estimates = c(a_norm = a_norm, sigma_mu = sigma_mu, b_norm = b_norm, sigma_pi = sigma_pi),
      posterior = data.frame(
        muS_posterior = muS_posterior,
        mu_posterior = mu_posterior,
        S_posterior = S_posterior,
        pi_posterior = pi_posterior
      )
    )
  return(ebe_normal)
}


# infer n (number of bins along one axis) from N_pairs = n*(n-1)/2
infer_n_from_long <- function(N_pairs) {
  n <- (-1 + sqrt(1 + 8 * N_pairs)) / 2
  n_int <- as.integer(round(n))
  if (n_int * (n_int + 1) / 2 != N_pairs)
    stop("N_pairs does not match n(n+1)/2; got N_pairs=", N_pairs, ", inferred n=", n_int)
  n_int
}

#' BZIP-GB(NB) imputation: muS (theta = mu*(1-S))
#' Also returns EB hyperparameters (a,b,c,d) per pair.
#'
#' @param sim_y       Integer matrix [N_pairs x K]
#' @param coverage    Numeric lambda multiplier (scalar OR length-K vector)
#' @param r           Neighborhood radius (default 1)
#' @param threshold   Borrow if (#nonzero in row) < threshold (default 3)
#' @param B           c + d (Beta concentration)
#' @param include_diag Logical, whether neighborhood includes diagonal (default TRUE)
#' @param return_hyper Logical, if TRUE returns list(muS, a, b, c, d)
#'
#' @return
#'   If return_hyper = FALSE (default):
#'       Numeric matrix muS [N_pairs x K]
#'   If return_hyper = TRUE:
#'       list(
#'         muS = [N_pairs x K] matrix,
#'         a   = numeric(N_pairs),
#'         b   = numeric(N_pairs),
#'         c   = numeric(N_pairs),
#'         d   = numeric(N_pairs)
#'       )
run_BZIP_GB_NB <- function(sim_y,
                           coverage,
                           r = 1,
                           threshold = 3,
                           B = 10,
                           include_diag = TRUE,
                           return_hyper = FALSE) {
  N <- nrow(sim_y); K <- ncol(sim_y)
  n_bins <- infer_n_from_long(N)
  neigh  <- get_neighbor_id(n_bins, r = r, include_diag = include_diag)
  
  muS <- matrix(0, nrow = N, ncol = K)
  colnames(muS) <- paste0("cell", seq_len(K))
  
  # storage for EB hyperparameters per pair
  a_vec <- numeric(N)
  b_mat <- matrix(NA, N, K)
  c_vec <- numeric(N)
  d_vec <- numeric(N)
  
  for (i in seq_len(N)) {
    Y_sim <- sim_y[i, ]
    
    # neighborhood borrowing
    if (sum(Y_sim > 0) >= threshold) {
      Y_input <- Y_sim
    } else {
      Y_input <- as.vector(sim_y[neigh[[i]], , drop = FALSE])
    }
    
    ebe <- get_EBE_ZNB_Gamma_Beta(Y_input, c_plus_d = B, lambda = coverage)
    
    # fallback if a is invalid
    if (!is.finite(ebe$a) || ebe$a < 0) {
      ebe_all <- get_EBE_ZNB_Gamma_Beta(as.vector(sim_y),
                                        c_plus_d = B,
                                        lambda   = coverage)
      ebe$a <- ebe_all$a
      ebe$b <- ebe_all$b
      # you can decide whether to also override c,d or keep pair-specific ones
      # here we keep original ebe$c, ebe$d unless they are missing
      if (is.null(ebe$c) || !is.finite(ebe$c) || ebe$c < 0) ebe$c <- ebe_all$c
      if (is.null(ebe$d) || !is.finite(ebe$d) || ebe$d < 0) ebe$d <- ebe_all$d
    }
    
    # store EB hyperparameters for this pair
    a_vec[i] <- ebe$a
    b_mat[i,] <- ebe$b
    c_vec[i] <- ebe$c
    d_vec[i] <- ebe$d
    
    # posterior pieces
    # E_S is vector length K if coverage is length-K
    E_S     <- get_E_S(ebe$a, ebe$b, ebe$c, ebe$d, coverage)  # E[S | Y=0]
    S_post  <- ifelse(Y_sim > 0, 0, E_S)
    mu_post <- (ebe$a + Y_sim) / (ebe$b + coverage)
    
    # theta = mu * (1 - S)
    muS[i, ] <- ifelse(Y_sim == 0, mu_post * (1 - E_S), mu_post)
  }
  
  if (return_hyper) {
    return(list(
      muS = muS,
      a   = a_vec,
      b   = b_mat,
      c   = c_vec,
      d   = d_vec
    ))
  } else {
    return(muS)
  }
}


#' Run BZIP-GB(NB) across multiple coverages
#' @param sim_y_list      list of integer matrices [N_pairs x K]
#' @param coverage_levels character OR numeric vector of coverages; same length/order as sim_y_list
#' @param r,threshold     neighborhood params
#' @param verbose         print progress and timing
#' @return                named list: names(coverage) -> muS matrix [N_pairs x K]
run_BZIP_GB_NB_list <- function(sim_y_list, coverage_levels, r = 1, threshold = 3, B = 10, include_diag = T, verbose = TRUE) {
  stopifnot(length(sim_y_list) == length(coverage_levels))
  cov_names <- as.character(coverage_levels)
  out <- setNames(vector("list", length(sim_y_list)), cov_names)
  
  for (i in seq_along(sim_y_list)) {
    cov_i <- coverage_levels[i]
    sim_y <- sim_y_list[[i]]
    
    if (verbose) cat(sprintf("Coverage = %s ... ", as.character(cov_i)))
    t0 <- Sys.time()
    muS_i <- run_BZIP_GB_NB(sim_y, coverage = as.numeric(cov_i), r = r, threshold = threshold, B = B, include_diag = include_diag)
    if (verbose) cat(sprintf("done in %.1f sec\n", as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    
    out[[i]] <- muS_i
  }
  out
}
