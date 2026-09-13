# Collect Higashi simulation exports into one manuscript-format RData object.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "tidyr", "jsonlite"))
source(file.path(script_dir, "simulation_benchmark_method_helpers.R"))

runs_root <- Sys.getenv("HICBZIP_SIM_HIGASHI_RUNS_ROOT", unset = path_here("results", "simulation", "benchmark_methods", "inputs", "Higashi", "runs"))
out_dir <- Sys.getenv("HICBZIP_SIM_HIGASHI_OUT_DIR", unset = path_here("results", "simulation", "benchmark_methods", "Higashi"))
ensure_dir(out_dir)

read_one_run <- function(run_dir) {
  metadata <- jsonlite::fromJSON(file.path(run_dir, "run_metadata.json"))
  export_dir <- file.path(run_dir, "export_long")
  files <- list(
    `Higashi(ori)` = list.files(export_dir, pattern = "^ori_.*\\.csv$", full.names = TRUE)[1],
    `Higashi(nbr0)` = list.files(export_dir, pattern = "^nbr0_.*\\.csv$", full.names = TRUE)[1],
    `Higashi(nbr5)` = list.files(export_dir, pattern = "^nbr5_.*\\.csv$", full.names = TRUE)[1]
  )
  if (any(is.na(unlist(files)))) stop("Missing Higashi exports in: ", export_dir, call. = FALSE)
  dplyr::bind_rows(lapply(names(files), function(model) {
    mat <- simulation_read_csv_matrix(files[[model]])
    tibble::tibble(model = model, chr = as.integer(metadata$chr), coverage = as.character(metadata$coverage),
                   muS = list(mat), N = nrow(mat), K = ncol(mat))
  }))
}

run_dirs <- list.dirs(runs_root, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^chr[0-9XY]+_cov", basename(run_dirs))]
if (length(run_dirs) == 0) stop("No Higashi simulation run directories found in: ", runs_root, call. = FALSE)

higashi_unified <- dplyr::bind_rows(lapply(sort(run_dirs), read_one_run)) |>
  dplyr::mutate(coverage = factor(.data$coverage, levels = simulation_coverages(), ordered = TRUE)) |>
  dplyr::arrange(.data$chr, .data$coverage, .data$model)

save(higashi_unified, file = file.path(out_dir, "higashi_unified_allchr.RData"))
message("Saved Higashi simulation object under: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
