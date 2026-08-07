
# res = My_HiCImpute(scHiC = single, 
#                    start_value = c(100,100,10,8,10,0.1,900,0.2,0,replicate(ncol(single),8)), 
#                    parameter = c(0.5, 5, 1.5), niter=1000, burnin=500, thin=2)



#########################################
#########################################

# mat_to_heatmap = function(mat, count_max = -1){
#   # Plot a heatmap, given a contact matrix. The entry value is bounded by count_max for visualization.
#   # mat_bounded = pmin(mat, count_max)
#   # heatmap(mat_bounded, Colv = NA, Rowv = NA, scale="none")
#   mat = as.matrix(mat)
#   mat[is.na(mat)] = 0
#   mat_long = expand.grid(bin1 = 1:nrow(mat), bin2 = 1:ncol(mat))
#   if (count_max < 0) count_max = max(mat)
#   mat_long$count = pmin(as.vector(mat), count_max)
# 
#   plt <- ggplot(mat_long, aes(bin1, bin2, fill=count))
#   plt <- plt + geom_tile()
#   plt <- plt + theme_minimal()
#   plt <- plt + scale_fill_gradient(low="white", high="blue", limits=c(0, count_max))
#   plt <- plt + xlab("") + ylab("")
#   plt
# }


simulate_mu_from_3Dcoords = function(coords,
                                     alpha0 = 5.6, alpha1 = -1,
                                     beta_l = 0.9, beta_g = 0.9, beta_m = 0.9){
  # Number of loci
  n <- nrow(coords)
  
  # Step 2: Generate the Euclidean 3D distance matrix
  distance_matrix <- as.matrix(dist(coords))
  
  # Step 3: Generate additional covariates (fragment length, GC content, and mappability)
  # set.seed(123456) # For reproducibility
  fragment_length <- runif(n, min=0.2, max=0.3) # Example range for fragment length
  gc_content <- runif(n, min=0.4, max=0.5) # GC content as a fraction
  mappability <- runif(n, min=0.9, max=1) # Mappability score
  
  # Coefficients
  alpha0 <- 5.6
  alpha1 <- -1
  beta_l <- beta_g <- beta_m <- 0.9
  
  # Using outer to compute the pairwise products for the covariates
  log_fragment_length_product <- outer(fragment_length, fragment_length, function(x, y) log(x * y))
  log_gc_content_product <- outer(gc_content, gc_content, function(x, y) log(x * y))
  log_mappability_product <- outer(mappability, mappability, function(x, y) log(x * y))
  
  # Using outer to compute the log of the distance matrix
  log_distance_matrix <- log(distance_matrix)
  
  # Compute the log_mu_matrix using vectorized operations
  log_mu_matrix <- alpha0 +
    alpha1 * log_distance_matrix +
    beta_l * log_fragment_length_product +
    beta_g * log_gc_content_product +
    beta_m * log_mappability_product
  
  # Exponentiate to get the mu_matrix
  mu_matrix <- exp(log_mu_matrix)
  
  # Set diagonal elements to 0 (assuming no self-interaction)
  sim_mu = mu_matrix[lower.tri(mu_matrix)]
  return(sim_mu)
}


mat_to_heatmap = function(mat, count_max = -1, bin_range = NULL){
  # Plot a heatmap, given a contact matrix. The entry value is bounded by count_max for visualization.
  # mat_bounded = pmin(mat, count_max)
  # heatmap(mat_bounded, Colv = NA, Rowv = NA, scale="none")
  mat = as.matrix(mat)
  mat[is.na(mat)] = 0
  if (!is.null(bin_range)) mat = mat[bin_range, bin_range]
  mat_long = expand.grid(bin1 = 1:nrow(mat), bin2 = 1:ncol(mat))
  if (count_max < 0) count_max = max(mat)
  mat_long$count = pmin(as.vector(mat), count_max)
  
  plt <- ggplot(mat_long, aes(bin1, bin2, fill=count))
  plt <- plt + geom_tile()
  plt <- plt + theme_minimal()
  plt <- plt + scale_fill_viridis_c(limits = c(0, count_max), option = "plasma")
  plt <- plt + xlab("") + ylab("") + theme(legend.position = "none")
  plt
}


mat_long_to_heatmaps = function(mat_long, n_plot, ncol = 3, maxi = -1, pi_S = F, bin_range = NULL,
                                off_diag = T, title = NULL, normalize = function(x) x, tick_by = 20){
  
  plots = apply(mat_long[,1:n_plot], 2, function(x){
    if (off_diag) x = matrix_long_to_matrix2D_offdiag(x) else x = matrix_long_to_matrix2D(x)
    x = normalize(x)
    p = mat_to_heatmap(x,maxi, bin_range)
    p = p + theme_minimal() + theme(legend.position = "none") + 
      scale_y_continuous(breaks=seq(0, nrow(x),by = tick_by)) + 
      scale_x_continuous(breaks=seq(0,nrow(x),by = tick_by))
    if (pi_S) p = mat_to_heatmap(x,bin_range = bin_range) + theme_minimal() + theme(legend.position = "none") + 
      scale_y_continuous(breaks=seq(0, nrow(x),by = tick_by)) + 
      scale_x_continuous(breaks=seq(0,nrow(x),by = tick_by)) + 
      scale_fill_gradient(low="white", high="red", name = "Probability", limits = c(0,1))
    return(p)
  })
  gridExtra::grid.arrange(grobs = plots, ncol = ncol, top = title)
}


mat_long_to_heatmap_compare = function(mat1, mat2){
  mat1 = matrix_long_to_matrix2D_offdiag(mat1)
  mat1 = log2(mat1+1)
  mat1 = mat1/max(mat1)
  
  mat2 = matrix_long_to_matrix2D_offdiag(mat2)
  mat2 = log2(mat2+1)
  mat2 = mat2/max(mat2)
  
  mat_plot = mat1
  mat_plot[lower.tri(mat_plot)] = mat2[lower.tri(mat2)]
  
  plt <- mat_to_heatmap(mat_plot) 
  return(plt)
}


pair_to_mtx = function(pair_file, size, res){
  # Given a pair library with two columns, output a 2D contact matrix
  n_bead_hic = ceiling(size/res)
  # n_bead_hic
  contact_mat = matrix(NA, n_bead_hic, n_bead_hic)
  bin_pair = ceiling(pair_file/res)
  colnames(bin_pair) = c("bin1", "bin2")
  
  bin_pair = bin_pair %>%
    mutate(count = 1) %>%
    group_by(bin1, bin2) %>%
    summarise(n = sum(count))
  
  for (i in 1:nrow(bin_pair)){
    contact_mat[bin_pair$bin1[i], bin_pair$bin2[i]] = bin_pair$n[i]
    contact_mat[bin_pair$bin2[i], bin_pair$bin1[i]] = bin_pair$n[i]
  }
  
  return(contact_mat)
}


matrix2D_to_matrix_long = function(mtx, include.diag = T){
  mtx_long = mtx[lower.tri(mtx, diag = include.diag)]
  mtx_long[is.na(mtx_long)] = 0
  return(mtx_long)
}


matrix3D_to_matrix_long = function(matrix3D){
  K = dim(matrix3D)[3]
  matrix_long = sapply(1:K, function(k){
    mtx = matrix3D[,,k]
    mtx_long = mtx[lower.tri(mtx, diag = T)]
    mtx_long[is.na(mtx_long)] = 0
    return(mtx_long)
  })
  return(matrix_long)
}


matrix_long_to_matrix2D = function(mtx_long){
  # Function to fill a square matrix with lower triangular values
  filllower <- function(lower_vals) {
    n <- (-1 + sqrt(1 + 8*length(lower_vals))) / 2
    if (floor(n) != n) {
      stop("Invalid number of lower triangular values!")
    }
    mat <- matrix(0, n, n)
    mat[lower.tri(mat, diag = T)] <- lower_vals
    return(mat)
  }
  
  # Fill the lower triangle
  lower_mat <- filllower(mtx_long)
  
  # Get full matrix
  full_mat <- lower_mat + t(lower_mat) - diag(diag(lower_mat))
  return(full_mat)
}

matrix_long_to_matrix2D = function(mtx_long, triangle = "lower"){
  # Function to fill a square matrix with lower triangular values
  filllower <- function(lower_vals) {
    n <- (-1 + sqrt(1 + 8*length(lower_vals))) / 2
    if (floor(n) != n) {
      stop("Invalid number of lower triangular values!")
    }
    mat <- matrix(0, n, n)
    if (triangle == "upper") mat[upper.tri(mat, diag = T)] <- lower_vals
    else mat[lower.tri(mat, diag = T)] <- lower_vals
    
    return(mat)
  }
  
  # Fill the lower triangle
  lower_mat <- filllower(mtx_long)
  
  # Get full matrix
  full_mat <- lower_mat + t(lower_mat)
  return(full_mat)
}


matrix_long_to_matrix2D_offdiag = function(mtx_long, triangle = "lower"){
  # Function to fill a square matrix with lower triangular values
  filllower <- function(lower_vals) {
    n <- (1 + sqrt(1 + 8*length(lower_vals))) / 2
    if (floor(n) != n) {
      stop("Invalid number of lower triangular values!")
    }
    mat <- matrix(0, n, n)
    if (triangle == "upper") mat[upper.tri(mat, diag = F)] <- lower_vals
    else mat[lower.tri(mat, diag = F)] <- lower_vals

    return(mat)
  }
  
  # Fill the lower triangle
  lower_mat <- filllower(mtx_long)
  
  # Get full matrix
  full_mat <- lower_mat + t(lower_mat)
  return(full_mat)
}

matrix_long_maxdist_to_matrix_long = function(mtx_long_md, n, idx_mtx_long, triangle = "lower"){
  mtx_long = rep(0, n*(n-1)/2)
  mtx_long[idx_mtx_long] = mtx_long_md
  return(mtx_long)
}

sim_scHiC_from_pair = function(pair_file, size, res, sparsity = 0.1){
  # pair_file = contact_pair_ch1; size = max(contact_pair_ch1); res = 16000; sparsity = 0.1;
  # set.seed(12345)
  index = sample(1:nrow(pair_file), ceiling(nrow(pair_file)*sparsity))
  pair_file_sm = pair_file[index]
  sc_mtx_sm = pair_to_mtx(pair_file_sm, size, res)
  return(sc_mtx_sm)
}


evaluate_mu_noscc = function(truth, mu){
  mse = mean((truth-mu)^2)
  mse_std = mean((truth-mu)^2) / mu^2
  return(c(mse = mse, mse_std = mse_std))
}

evaluate_mu = function(truth, mu, smooth = 5, res = 50000){
  cor = cor(truth, mu)
  cor_log = cor(log(truth+1), log(mu+1))
  mse = mean((truth-mu)^2)
  mse_std = mean((truth-mu)^2) / var(truth)
  m1 = matrix_long_to_matrix2D_offdiag(truth)
  m2 = matrix_long_to_matrix2D_offdiag(mu)
  scc.out = hicrep::get.scc(m1, m2,
                            res, smooth, 0, min(5e+06, res*(nrow(m1)-1)))
  return(c(cor = cor, cor_log = cor_log,
           mse = mse, mse_std = mse_std,
           scc = scc.out$scc, scc_sd = scc.out$std))
}


check_heatmap = function(impute_mu, impute_pi, k,
                         single, truth, cell_name, 
                         pi_same_count, max_count = 10){
  
  p_truth = truth[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    .[bin_range, bin_range] %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("True Intensity of cell #", k, " from ", cell_name, " SCL Structure"), subtitle = plot_subtitle)
  
  p_single = single[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Simulated scHiC of cell #", k, " from ", cell_name, " SCL Structure"), subtitle = plot_subtitle)
  
  p_pi_0.5 = (impute_mu[,k] * (impute_pi < 0.5)) %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Imputed Contact Intensity of cell #", k, " from ", cell_name," (pi < 0.5)"), subtitle = plot_subtitle)
  
  p_pi_mean = (impute_mu[,k] * (impute_pi < mean(impute_pi))) %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., 10) +
    labs(title = paste0("Imputed Contact Intensity of cell #", k, " from ", cell_name," (pi < mean)"), subtitle = plot_subtitle)
  
  p_pi_samecount = (impute_mu[,k] * (impute_pi < pi_same_count[k])) %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., 10) +
    labs(title = paste0("Imputed Contact Intensity of cell #", k, " from ", cell_name," (pi s.t. same count)"), subtitle = plot_subtitle)
  
  
  p_mu = impute_mu[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Posterior mu of cell #", k," from ", cell_name), subtitle = plot_subtitle)
  
  p_pi = impute_pi %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Posterior of pi of cell #", k," from ", cell_name), subtitle = plot_subtitle) +
    scale_fill_gradient(low="white", high="red", name = "Probability", limits = c(0,1))
  
  p_mu_pi = (impute_mu[,k]*(1-impute_pi)) %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., 10) +
    labs(title = paste0("Posterior mu*(1-pi) of cell #", k," from ", cell_name), subtitle = plot_subtitle)
  
  return(p_list = list(p_truth, p_single, p_mu, p_pi, p_mu_pi, p_pi_0.5, p_pi_mean, p_pi_samecount))
}



check_heatmap_sim_mu_pi = function(impute_mu, impute_pi, k, single,
                                   CTCF_prior, true_mu, true_pi, cell_name, 
                                   max_count = 10){
  
  p_truth_mu = true_mu[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    .[bin_range, bin_range] %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("True Mu of cell #", k, " from ", cell_name, " SCL Structure"), subtitle = plot_subtitle)
  
  p_truth_pi = true_pi %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    .[bin_range, bin_range] %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("True Pi of cell #", k, " from ", cell_name, " SCL Structure"), subtitle = plot_subtitle) +
    scale_fill_gradient(low="white", high="red", name = "Probability", limits = c(0,1))
  
  p_single = single[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Simulated scHiC of cell #", k, " from ", cell_name, " SCL Structure"), subtitle = plot_subtitle)
  
  # p_pi_0.5 = (impute_mu[,k] * (impute_pi < 0.5)) %>%
  #   matrix_long_to_matrix2D_offdiag(.) %>%
  #   mat_to_heatmap(., max_count) +
  #   labs(title = paste0("Imputed Contact Intensity of cell #", k, " from ", cell_name," (pi < 0.5)"), subtitle = plot_subtitle)
  # 
  # p_pi_mean = (impute_mu[,k] * (impute_pi < mean(impute_pi))) %>%
  #   matrix_long_to_matrix2D_offdiag(.) %>%
  #   mat_to_heatmap(., 10) +
  #   labs(title = paste0("Imputed Contact Intensity of cell #", k, " from ", cell_name," (pi < mean)"), subtitle = plot_subtitle)
  # 
  # p_pi_samecount = (impute_mu[,k] * (impute_pi < pi_same_count[k])) %>%
  #   matrix_long_to_matrix2D_offdiag(.) %>%
  #   mat_to_heatmap(., 10) +
  #   labs(title = paste0("Imputed Contact Intensity of cell #", k, " from ", cell_name," (pi s.t. same count)"), subtitle = plot_subtitle)
  
  p_prior = CTCF_prior %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(.) +
    labs(title = "CTCF Prior on mu (win.size=5M)", subtitle = plot_subtitle)
    
  p_mu = impute_mu[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Posterior mu of cell #", k," from ", cell_name), subtitle = plot_subtitle)
  
  p_pi = impute_pi %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Posterior of pi of cell #", k," from ", cell_name), subtitle = plot_subtitle) +
    scale_fill_gradient(low="white", high="red", name = "Probability", limits = c(0,1))
  
  # p_mu_pi = (impute_mu[,k]*(1-impute_pi)) %>%
  #   matrix_long_to_matrix2D_offdiag(.) %>%
  #   mat_to_heatmap(., 10) +
  #   labs(title = paste0("Posterior mu*(1-pi) of cell #", k," from ", cell_name), subtitle = plot_subtitle)
  
  return(p_list = list(p_single, p_truth_mu, p_truth_pi, p_prior, p_mu, p_pi))
}



check_heatmap_sim_mu_pi_y_s = function(impute_mu, impute_pi, k, single,
                                   CTCF_prior, true_mu, true_pi, true_s, true_y,
                                   cell_name, pi_same_count, max_count = 10){
  
  p_truth_mu = true_mu[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    .[bin_range, bin_range] %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("True Mu of cell #", k, " from ", cell_name, " SCL Structure"), subtitle = plot_subtitle)
  
  p_truth_pi = true_pi %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    .[bin_range, bin_range] %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("True Pi of cell #", k, " from ", cell_name, " SCL Structure"), subtitle = plot_subtitle) +
    scale_fill_gradient(low="white", high="red", name = "Probability", limits = c(0,1))
  
  p_truth_s = true_s[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    .[bin_range, bin_range] %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("True SZ of cell #", k, " from ", cell_name, " SCL Structure"), subtitle = plot_subtitle) +
    scale_fill_gradient(low="white", high="red", name = "Probability", limits = c(0,1))
  
  p_truth_y = true_y[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    .[bin_range, bin_range] %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("True Contact Intensity of cell #", k, " from ", cell_name, " SCL Structure"), subtitle = plot_subtitle)
  
  p_single = single[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Simulated scHiC of cell #", k, " from ", cell_name, " SCL Structure"), subtitle = plot_subtitle)
  
  p_pi_0.5 = (impute_mu[,k] * (impute_pi < 0.5)) %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Imputed Contact Intensity of cell #", k, " from ", cell_name," (pi < 0.5)"), subtitle = plot_subtitle)

  p_pi_mean = (impute_mu[,k] * (impute_pi < mean(impute_pi))) %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., 10) +
    labs(title = paste0("Imputed Contact Intensity of cell #", k, " from ", cell_name," (pi < mean)"), subtitle = plot_subtitle)

  p_pi_samecount = (impute_mu[,k] * (impute_pi < pi_same_count[k])) %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., 10) +
    labs(title = paste0("Imputed Contact Intensity of cell #", k, " from ", cell_name," (pi s.t. same count)"), subtitle = plot_subtitle)
  
  p_prior = CTCF_prior %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(.) +
    labs(title = "CTCF Prior on mu (win.size=5M)", subtitle = plot_subtitle)
  
  p_mu = impute_mu[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Posterior mu of cell #", k," from ", cell_name), subtitle = plot_subtitle)
  
  p_pi = impute_pi[,k] %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., max_count) +
    labs(title = paste0("Posterior of pi of cell #", k," from ", cell_name), subtitle = plot_subtitle) +
    scale_fill_gradient(low="white", high="red", name = "Probability", limits = c(0,1))
  
  p_mu_pi = (impute_mu[,k]*(1-impute_pi)) %>%
    matrix_long_to_matrix2D_offdiag(.) %>%
    mat_to_heatmap(., 10) +
    labs(title = paste0("Posterior mu*(1-pi) of cell #", k," from ", cell_name), subtitle = plot_subtitle)
  
  return(p_list = list(p_single, p_truth_mu, p_truth_pi, p_truth_y,
                       p_prior, p_mu, p_pi, p_pi_samecount))
}


fit_TopDom = function(mat_long, start, end, resolution){
  muS_mtx_2D = data.frame(
    chr = "chr1",
    from.coord = seq(start, end-resolution, by = resolution),
    to.coord = seq(start+resolution, end, by = resolution)
  )
  muS_mtx_2D = cbind(muS_mtx_2D, matrix_long_to_matrix2D_offdiag(mat_long))
  write.table(muS_mtx_2D, file = "temp.txt", sep = " ", row.names = F, col.names = F)
  data_true_muS <- readHiC("temp.txt")
  fit <- TopDom(data_true_muS, window.size = 5L)
  return(fit)
}

Plot_TopDom_with_Truth = function(mat_long, truth, start, end, resolution,
                                  normalize = function(x) log2(x+1)/log2(max(x)+1), title = NULL){
  fit_true <- fit_TopDom(truth, start, end, resolution)
  fit <- fit_TopDom(mat_long, start, end, resolution)
  
  
  mat1 = matrix_long_to_matrix2D_offdiag(truth)
  mat1 = normalize(mat1)
  
  mat2 = matrix_long_to_matrix2D_offdiag(mat_long)
  mat2 = normalize(mat2)
  
  mat_plot = mat1
  mat_plot[lower.tri(mat_plot)] = mat2[lower.tri(mat2)]
  
  pp <- mat_to_heatmap(mat_plot) 
  for (kk in seq_len(nrow(fit_true$domain))) {
    x0 <- fit_true$domain$from.id[kk]
    y0 <- fit_true$domain$to.id[kk]
    pp <- pp + geom_segment(x = x0, y = y0, xend = x0, yend = x0, color = "green", linewidth=1) +
      geom_segment(x = x0, y = y0, xend = y0, yend = y0, color = "green", linewidth=1)
  }
  for (kk in seq_len(nrow(fit$domain))) {
    y0 <- fit$domain$from.id[kk]
    x0 <- fit$domain$to.id[kk]
    pp <- pp + geom_segment(x = x0, y = y0, xend = x0, yend = x0, color = "red", linewidth=1) +
      geom_segment(x = x0, y = y0, xend = y0, yend = y0, color = "red", linewidth=1)
  }
  pp <- pp + theme(legend.position = "none") + ggtitle(title)
  return(pp)
}


simulate_mu_from_3Dcoords = function(coords,
                                     alpha0 = 5.6, alpha1 = -1,
                                     beta_l = 0.9, beta_g = 0.9, beta_m = 0.9){
  # Number of loci
  n <- nrow(coords)
  
  # Step 2: Generate the Euclidean 3D distance matrix
  distance_matrix <- as.matrix(dist(coords))
  
  # Step 3: Generate additional covariates (fragment length, GC content, and mappability)
  # set.seed(123456) # For reproducibility
  fragment_length <- runif(n, min=0.2, max=0.3) # Example range for fragment length
  gc_content <- runif(n, min=0.4, max=0.5) # GC content as a fraction
  mappability <- runif(n, min=0.9, max=1) # Mappability score
  
  # Coefficients
  alpha0 <- 5.6
  alpha1 <- -1
  beta_l <- beta_g <- beta_m <- 0.9
  
  # Using outer to compute the pairwise products for the covariates
  log_fragment_length_product <- outer(fragment_length, fragment_length, function(x, y) log(x * y))
  log_gc_content_product <- outer(gc_content, gc_content, function(x, y) log(x * y))
  log_mappability_product <- outer(mappability, mappability, function(x, y) log(x * y))
  
  # Using outer to compute the log of the distance matrix
  log_distance_matrix <- log(distance_matrix)
  
  # Compute the log_mu_matrix using vectorized operations
  log_mu_matrix <- alpha0 +
    alpha1 * log_distance_matrix +
    beta_l * log_fragment_length_product +
    beta_g * log_gc_content_product +
    beta_m * log_mappability_product
  
  # Exponentiate to get the mu_matrix
  mu_matrix <- exp(log_mu_matrix)
  
  # Set diagonal elements to 0 (assuming no self-interaction)
  sim_mu = mu_matrix[lower.tri(mu_matrix)]
  return(sim_mu)
}

process_hdf5_file <- function(impute_file, include_diag = F) {
  # Open the HDF5 file
  h5f <- H5Fopen(impute_file)
  
  # Extract the CSR components
  data <- h5f$Matrix$data %>% as.vector()
  indices <- h5f$Matrix$indices %>% as.vector()
  indptr <- h5f$Matrix$indptr %>% as.vector()
  H5Fclose(h5f)
  
  # Determine the number of rows and columns
  num_rows <- length(indptr) - 1
  num_cols <- max(indices) + 1
  
  # Create the sparse matrix in CSR format
  sparse_matrix <- sparseMatrix(
    i = rep(1:num_rows, diff(indptr)),
    j = indices + 1,  # Convert from 0-based to 1-based indexing
    x = data,
    dims = c(num_rows, num_cols)
  )
  
  # Convert the sparse matrix to a dense matrix
  dense_matrix <- as.matrix(sparse_matrix)
  dense_matrix <- dense_matrix[,-1]
  dense_matrix <- dense_matrix[-1,]
  dense_matrix[lower.tri(dense_matrix, diag = include_diag)] = t(dense_matrix)[lower.tri(t(dense_matrix), diag = include_diag)]
  
  # Convert the dense matrix to long format
  dense_matrix_long <- as.vector(dense_matrix[lower.tri(dense_matrix, diag = include_diag)])
  # dense_matrix_long %>% matrix_long_to_matrix2D_offdiag() %>% mat_to_heatmap()
  
  return(dense_matrix_long)
}

draw_traceplot = function(impute_full, par, iter_id, truth = NULL, coverage, save = F, plot_par = NULL){
  stan_df <- as.data.frame(impute_full[,,par]) %>%
    mutate(iteration = row_number()) %>%     # Add iteration number
    pivot_longer(cols = -iteration, names_to = "Method", values_to = "value") %>%
    group_by(Method) %>%
    mutate(cum_mean = cummean(value),
           Method = str_replace(Method, "chain", "Chain (STAN)"))
  if (is.null(plot_par)) plot_par = par
  p1 = ggplot(rbind(stan_df %>% select(-cum_mean)), aes(x = iteration, y = value, color = Method)) +
    geom_line(alpha = 0.5) +
    labs(x = "Iteration", y = "Value", title = paste("Trace Plot of", plot_par)) +
    theme_minimal() 
  stan_df <- as.data.frame(impute_full[iter_id,,par]) %>%
    mutate(iteration = row_number()) %>%     # Add iteration number
    pivot_longer(cols = -iteration, names_to = "Method", values_to = "value") %>%
    group_by(Method) %>%
    mutate(cum_mean = cummean(value),
           Method = str_replace(Method, "chain", "Chain (STAN)"))
  p2 = ggplot(rbind(stan_df %>% select(-value)), aes(x = iteration, y = cum_mean, color = Method)) +
    geom_line() +
    labs(x = "Iteration", y = "Cumulative Mean", title =  paste("Cumulative Mean of", plot_par)) +
    theme_minimal() 
  if (!is.null(truth)){
    p1 = p1 + geom_hline(yintercept = truth*coverage, linetype = "dashed", color = "red")
    p2 = p2 + geom_hline(yintercept = truth*coverage, linetype = "dashed", color = "red")
  }
  if (save){
    png(paste0("tr_", gsub("[\\[\\],]", "_", par, perl = TRUE), ".png"), 1600, 200)
    grid.arrange(p1,p2,ncol = 2)
    dev.off()
  }else{
    grid.arrange(p1,p2,ncol =1)
  }
}


# Function to calculate insulation scores for Hi-C contact matrix
calculate_insulation_score <- function(contact_matrix, window_size = 5, bin_size = 50*10e3) {
  # Get matrix dimensions
  n_bins <- nrow(contact_matrix)
  
  # Initialize vectors to store results
  IS <- numeric(n_bins)
  IS_prime <- numeric(n_bins)
  
  # Calculate the number of bins to exclude at the start and end
  n_exclude <- window_size
  
  # Calculate raw insulation score (IS) for each bin
  for(i in (n_exclude + 1):(n_bins - n_exclude)) {
    # Define the bin sets J and K
    J_bins <- (i - window_size):(i - 1)  # [i-5, i-1]
    K_bins <- (i + 1):(i + window_size)  # [i+1, i+5]
    
    # Sum the contacts between bins in J and K
    contact_sum <- 0
    for(j in J_bins) {
      for(k in K_bins) {
        contact_sum <- contact_sum + contact_matrix[j, k]
      }
    }
    
    # Calculate IS for this bin
    # n^2 is window_size^2 since n = window_size
    IS[i] <- contact_sum / (window_size * window_size)
  }
  
  # Calculate the average IS (ISavg)
  valid_IS <- IS[IS != 0]  # Exclude zeros from boundary regions
  IS_avg <- mean(valid_IS, na.rm = TRUE)
  
  # Calculate IS' using log2 transformation
  for(i in 1:n_bins) {
    if(IS[i] != 0) {
      IS_prime[i] <- log2((IS[i] / IS_avg) + 1)
    } else {
      IS_prime[i] <- NA  # Set boundary regions to NA
    }
  }
  
  # Create a data frame with results
  results <- data.frame(
    bin_index = 1:n_bins,
    position = (1:n_bins) * bin_size,
    IS = IS,
    IS_prime = IS_prime
  )
  
  return(results)
}

get_IS_cor = function(true, impute, n_cell, win_size, resolution){
  sapply(1:n_cell, function(k){
    IS_true = true[,k] %>% matrix_long_to_matrix2D_offdiag() %>% calculate_insulation_score(., window_size = win_size, bin_size = resolution)
    IS_impute = impute[,k] %>% matrix_long_to_matrix2D_offdiag() %>% calculate_insulation_score(., window_size = win_size, bin_size = resolution)
    if (sum(!is.na(IS_impute$IS_prime)) == 0) return(NA)
    else return(cor(IS_true$IS_prime, IS_impute$IS_prime, use="complete.obs"))
  })
}


get_IS_mse <- function(true, impute, n_cell, win_size, resolution) {
  out <- lapply(seq_len(n_cell), function(k) {
    
    IS_true   <- true[, k]   %>% matrix_long_to_matrix2D_offdiag() %>%
      calculate_insulation_score(window_size = win_size, bin_size = resolution)
    
    IS_impute <- impute[, k] %>% matrix_long_to_matrix2D_offdiag() %>%
      calculate_insulation_score(window_size = win_size, bin_size = resolution)
    
    # align + keep complete cases
    ok <- is.finite(IS_true$IS_prime) & is.finite(IS_impute$IS_prime)
    
    if (sum(ok) == 0) {
      return(c(MSE = NA_real_, sMSE = NA_real_))
    }
    
    diff <- IS_impute$IS_prime[ok] - IS_true$IS_prime[ok]
    mse  <- mean(diff^2)
    
    vtrue <- stats::var(IS_true$IS_prime[ok])
    smse  <- if (is.na(vtrue) || vtrue == 0) NA_real_ else mse / vtrue
    
    c(MSE = mse, sMSE = smse)
  })
  
  do.call(rbind, out)  # n_cell x 2 matrix
}