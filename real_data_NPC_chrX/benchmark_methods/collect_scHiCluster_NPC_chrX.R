# Collect scHiCluster NPC chrX HDF5 outputs into the common long-vector format.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "stringr", "rhdf5", "Matrix"))
source_hicbzip_core()
source(file.path(script_dir, "..", "NPC_chrX_workflow_helpers.R"))

input_root <- Sys.getenv(
  "HICBZIP_REAL1_SCHICLUSTER_ROOT",
  unset = path_here("results", "NPC_chrX", "benchmark_inputs", "scHiCluster")
)
out_dir <- path_here("results", "NPC_chrX", "scHiCluster")
out_file <- file.path(out_dir, "o.hicluster.muS.full.AllCoverage.RData")

require_dirs(input_root, label = "scHiCluster NPC chrX output root")
ensure_dir(out_dir)

coverage_labels <- npc_chrX_coverage_labels()
o.hicluster.muS.full <- lapply(coverage_labels, function(cov_label) {
  cov_dir <- file.path(input_root, paste0("NPC250k_0h_X_full_", cov_label), "imputed_matrix")
  require_dirs(cov_dir, label = paste("scHiCluster output for coverage", cov_label))
  h5_files <- list.files(cov_dir, pattern = "_chrX\\.hdf5$", full.names = TRUE)
  if (length(h5_files) == 0) stop("No scHiCluster HDF5 files found in: ", cov_dir, call. = FALSE)
  h5_files <- sort(h5_files)
  long_list <- lapply(h5_files, process_hdf5_file, include_diag = TRUE)
  do.call(cbind, long_list)
})
names(o.hicluster.muS.full) <- paste0("scHiCluster,ld=", coverage_labels)

save(o.hicluster.muS.full, file = out_file)
message("Saved: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
