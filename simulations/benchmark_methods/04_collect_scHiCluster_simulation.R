# Collect scHiCluster simulation HDF5 outputs into one manuscript-format RDS object.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "rhdf5", "Matrix"))
source_hicbzip_core()
source(file.path(script_dir, "simulation_benchmark_method_helpers.R"))

input_root <- Sys.getenv(
  "HICBZIP_SIM_SCHICLUSTER_ROOT",
  unset = path_here("results", "simulation", "benchmark_methods", "inputs", "scHiCluster")
)
out_dir <- Sys.getenv(
  "HICBZIP_SIM_SCHICLUSTER_OUT_DIR",
  unset = path_here("results", "simulation", "benchmark_methods", "scHiCluster")
)
ensure_dir(out_dir)
coverage_labels <- simulation_coverages()

num_sort <- function(paths) {
  idx <- suppressWarnings(as.numeric(gsub(".*_([0-9]+)_chr[0-9XY]+\\.hdf5$", "\\1", basename(paths))))
  paths[order(idx, na.last = TRUE)]
}

schicluster_muS_all <- list()
for (coverage_label in coverage_labels) {
  imp_dir <- file.path(input_root, paste0("sim_HBA3_random_chr_", coverage_label), "imputed_matrix")
  require_dirs(imp_dir, label = paste("scHiCluster simulation output for coverage", coverage_label))
  h5_files <- list.files(imp_dir, pattern = "\\.hdf5$", full.names = TRUE)
  if (length(h5_files) == 0) stop("No scHiCluster HDF5 files found in: ", imp_dir, call. = FALSE)
  chr_tags <- unique(sub(".*_chr([0-9XY]+)\\.hdf5$", "chr\\1", basename(h5_files)))

  for (chr_tag in chr_tags) {
    chr_files <- h5_files[grepl(paste0("_", chr_tag, "\\.hdf5$"), h5_files)]
    chr_files <- num_sort(chr_files)
    long_list <- lapply(chr_files, process_hdf5_file, include_diag = FALSE)
    muS_mat <- do.call(cbind, long_list)
    colnames(muS_mat) <- paste0("cell", seq_along(long_list))
    if (is.null(schicluster_muS_all[[chr_tag]])) schicluster_muS_all[[chr_tag]] <- list()
    schicluster_muS_all[[chr_tag]][[coverage_label]] <- muS_mat
  }
}

saveRDS(schicluster_muS_all, file.path(out_dir, "scHiCluster_imputed_muS_allchr_allcov.rds"))
message("Saved scHiCluster simulation object under: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
