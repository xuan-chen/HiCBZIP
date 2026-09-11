# Collect Higashi NPC chrX exports into the common long-vector format.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "tidyr", "stringr"))
source(file.path(script_dir, "..", "NPC_chrX_workflow_helpers.R"))

runs_root <- Sys.getenv(
  "HICBZIP_REAL1_HIGASHI_RUNS_ROOT",
  unset = path_here("results", "NPC_chrX", "benchmark_inputs", "Higashi", "runs_chrX")
)
out_dir <- path_here("results", "NPC_chrX", "Higashi")
out_file <- file.path(out_dir, "higashi_unified_chrX_NPC250k_diag.RData")

require_dirs(runs_root, label = "Higashi NPC chrX runs root")
ensure_dir(out_dir)

read_one_run <- function(run_dir) {
  cov <- sub("^chrX_cov_", "", basename(run_dir))
  cov_label <- npc_chrX_normalize_cov_label(cov)
  export_dir <- file.path(run_dir, "export_long_diag")
  require_dirs(export_dir, label = paste("Higashi export for coverage", cov))
  files <- list(
    `Higashi(ori)` = list.files(export_dir, pattern = "^ori_diag_.*\\.csv$", full.names = TRUE)[1],
    `Higashi(nbr0)` = list.files(export_dir, pattern = "^nbr0_diag_.*\\.csv$", full.names = TRUE)[1],
    `Higashi(nbr5)` = list.files(export_dir, pattern = "^nbr5_diag_.*\\.csv$", full.names = TRUE)[1]
  )
  if (any(is.na(unlist(files)))) stop("Missing one or more Higashi CSV exports in: ", export_dir, call. = FALSE)
  dplyr::bind_rows(lapply(names(files), function(model) {
    mat <- npc_chrX_read_csv_matrix(files[[model]])
    tibble::tibble(model = model, chr = "X", coverage = cov_label, muS = list(mat), N = nrow(mat), K = ncol(mat))
  }))
}

run_dirs <- list.dirs(runs_root, recursive = FALSE, full.names = TRUE)
run_dirs <- run_dirs[grepl("^chrX_cov_", basename(run_dirs))]
if (length(run_dirs) == 0) stop("No Higashi run directories found in: ", runs_root, call. = FALSE)

higashi_unified_chrX <- dplyr::bind_rows(lapply(sort(run_dirs), read_one_run)) |>
  dplyr::mutate(coverage = factor(.data$coverage, levels = npc_chrX_coverage_labels(), ordered = TRUE)) |>
  dplyr::arrange(.data$chr, .data$coverage, .data$model)

save(higashi_unified_chrX, file = out_file)
message("Saved: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
