# Collect Fast-Higashi simulation exports into one manuscript-format RData object.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "tidyr", "jsonlite"))
source(file.path(script_dir, "simulation_benchmark_method_helpers.R"))

runs_root <- Sys.getenv("HICBZIP_SIM_FASTHIGASHI_RUNS_ROOT", unset = path_here("results", "simulation", "benchmark_methods", "inputs", "FastHigashi", "runs"))
out_dir <- Sys.getenv("HICBZIP_SIM_FASTHIGASHI_OUT_DIR", unset = path_here("results", "simulation", "benchmark_methods", "FastHigashi"))
ensure_dir(out_dir)

read_one_run <- function(run_dir) {
  metadata <- jsonlite::fromJSON(file.path(run_dir, "run_metadata.json"))
  export_dir <- file.path(run_dir, "export_long")
  csv_file <- list.files(export_dir, pattern = "^fasthigashi_prwr_.*\\.csv$", full.names = TRUE)[1]
  if (is.na(csv_file)) stop("Missing Fast-Higashi export in: ", export_dir, call. = FALSE)
  mat <- simulation_read_csv_matrix(csv_file)
  tibble::tibble(
    model = "Fast-Higashi(prwr,rwr_covlt0.1conv,rank3)",
    chr = as.integer(metadata$chr),
    coverage = as.character(metadata$coverage),
    muS = list(mat),
    N = nrow(mat),
    K = ncol(mat)
  )
}

run_dirs <- list.dirs(runs_root, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^chr[0-9XY]+_cov", basename(run_dirs))]
if (length(run_dirs) == 0) stop("No Fast-Higashi simulation run directories found in: ", runs_root, call. = FALSE)

fasthigashi_unified <- dplyr::bind_rows(lapply(sort(run_dirs), read_one_run)) |>
  dplyr::mutate(coverage = factor(.data$coverage, levels = simulation_coverages(), ordered = TRUE)) |>
  dplyr::arrange(.data$chr, .data$coverage)

save(fasthigashi_unified, file = file.path(out_dir, "fasthigashi_unified_allchr_rank3.RData"))
message("Saved Fast-Higashi simulation object under: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
