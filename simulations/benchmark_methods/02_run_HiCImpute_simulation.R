# Run HiCImpute for the manuscript simulation benchmark.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
require_packages(c("HiCImpute"))
source_hicbzip_core()
source(file.path(script_dir, "simulation_benchmark_method_helpers.R"))

input_dir <- Sys.getenv("HICBZIP_SIM_INPUT_DIR", unset = path_here("data", "processed", "simulation", "input"))
out_dir <- Sys.getenv("HICBZIP_SIM_HICIMPUTE_OUT_DIR", unset = path_here("results", "simulation", "benchmark_methods", "HiCImpute"))
ensure_dir(out_dir)

coverage_labels <- simulation_coverages()
sim_files <- simulation_list_input_files(input_dir)
if (length(sim_files) == 0) stop("No processed simulation input files found in: ", input_dir, call. = FALSE)

results_by_chr <- list()
for (input_file in sim_files) {
  env <- new.env(parent = emptyenv())
  load(input_file, envir = env)
  if (!exists("true_muS", envir = env)) stop("Missing true_muS in: ", input_file, call. = FALSE)
  true_muS <- as.matrix(env$true_muS)
  chr <- paste0("chr", simulation_extract_chr(input_file))
  n_bins <- simulation_infer_n_bins_offdiag(nrow(true_muS))
  results_by_chr[[chr]] <- list()

  for (coverage_label in coverage_labels) {
    set.seed(123456)
    sim_y <- matrix(
      rpois(nrow(true_muS) * ncol(true_muS), lambda = as.numeric(coverage_label) * as.numeric(true_muS)),
      nrow = nrow(true_muS),
      ncol = ncol(true_muS)
    )
    scHiC <- apply(sim_y, 2, function(y) {
      m2d <- matrix_long_to_matrix2D_offdiag(y)
      m2d[upper.tri(m2d, diag = FALSE)]
    })
    o.hicimpute <- HiCImpute::MCMCImpute(
      scHiC = scHiC,
      bulk = rowSums(scHiC),
      expected = NULL,
      startval = c(10, 1, 10, 1, 10, 0.5, 10, 0.5, 0, replicate(ncol(scHiC), 1)),
      n = n_bins,
      mc.cores = as.integer(Sys.getenv("HICBZIP_HICIMPUTE_WORKERS", unset = "1")),
      cutoff = 0.5,
      niter = as.integer(Sys.getenv("HICBZIP_HICIMPUTE_NITER", unset = "1000")),
      burnin = as.integer(Sys.getenv("HICBZIP_HICIMPUTE_BURNIN", unset = "500"))
    )
    muS_est <- apply(o.hicimpute$Impute_SZ, 2, function(x) {
      m2d <- matrix_long_to_matrix2D_offdiag(x, triangle = "upper")
      matrix2D_to_matrix_long(m2d, include.diag = FALSE)
    })
    results_by_chr[[chr]][[coverage_label]] <- muS_est
  }
  saveRDS(results_by_chr[[chr]], file.path(out_dir, paste0("results_", chr, ".rds")))
}

saveRDS(results_by_chr, file.path(out_dir, "HiCImpute_muS_allchr_allcov.rds"))
message("Saved HiCImpute simulation outputs under: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
