required_score_pkgs <- c("jsonlite", "dplyr", "tibble", "tidyr")
missing_score_pkgs <- required_score_pkgs[
  !vapply(required_score_pkgs, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_score_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_score_pkgs, collapse = ", "), call. = FALSE)
}

library(jsonlite)
library(dplyr)
library(tibble)
library(tidyr)

score_metric_keys <- c(
  "ari_k-means", "ari_agglomerative", "ari_gmm", "ari_louvain", "ari_leiden", "best_ari",
  "nmi_k-means", "nmi_agglomerative", "nmi_gmm", "nmi_louvain", "nmi_leiden", "best_nmi",
  "silhouette_k-means", "silhouette_agglomerative", "silhouette_gmm",
  "silhouette_louvain", "silhouette_leiden", "best_silhouette",
  "best_silhouette-gt", "wall_time"
)

resolve_score_bin <- function() {
  score_bin <- Sys.which("score")
  if (!nzchar(score_bin)) {
    candidates <- c(Sys.getenv("HICBZIP_SCORE_BIN", unset = ""))
    hit <- candidates[file.exists(candidates)]
    if (length(hit) > 0) score_bin <- hit[[1]]
  }
  if (!nzchar(score_bin)) {
    stop("`score` command not found. Activate the SCORE environment or set HICBZIP_SCORE_BIN.", call. = FALSE)
  }
  score_bin
}

run_score_cmd <- function(score_bin, args, label, stdout_file = NULL, stderr_file = NULL, force_cpu = FALSE) {
  cat("\n[", label, "]\n", sep = "")
  if (force_cpu) cat("CUDA_VISIBLE_DEVICES=\"\" ")
  cat("score ", paste(args, collapse = " "), "\n", sep = "")

  if (is.null(stdout_file)) stdout_file <- tempfile(fileext = ".stdout.log")
  if (is.null(stderr_file)) stderr_file <- tempfile(fileext = ".stderr.log")

  old_cuda <- Sys.getenv("CUDA_VISIBLE_DEVICES", unset = NA_character_)
  if (force_cpu) Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  on.exit({
    if (force_cpu) {
      if (is.na(old_cuda)) Sys.unsetenv("CUDA_VISIBLE_DEVICES") else Sys.setenv(CUDA_VISIBLE_DEVICES = old_cuda)
    }
  }, add = TRUE)

  status <- system2(score_bin, args = args, stdout = stdout_file, stderr = stderr_file)
  stdout_lines <- if (file.exists(stdout_file)) readLines(stdout_file, warn = FALSE) else character(0)
  stderr_lines <- if (file.exists(stderr_file)) readLines(stderr_file, warn = FALSE) else character(0)

  cat("Return code:", status, "\n")
  if (length(stdout_lines) > 0) cat(paste(stdout_lines, collapse = "\n"), "\n")
  if (!identical(status, 0L)) {
    cat("\n--- STDERR ---\n")
    if (length(stderr_lines) > 0) cat(paste(stderr_lines, collapse = "\n"), "\n")
    stop(label, " failed.", call. = FALSE)
  }

  invisible(list(status = as.integer(status), stdout = stdout_lines, stderr = stderr_lines))
}

extract_score_metric <- function(path, key) {
  if (!file.exists(path)) return(NA_real_)
  x <- jsonlite::fromJSON(path, simplifyVector = FALSE)[[key]]
  if (is.null(x)) return(NA_real_)
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  vals <- suppressWarnings(as.numeric(x))
  if (length(vals) == 0 || all(is.na(vals))) return(NA_real_)
  vals[which(!is.na(vals))[1]]
}

summarize_score_results <- function(results_df, group_cols = "method") {
  results_df %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(score_metric_keys),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::group_by(dplyr::across(tidyselect::all_of(c(group_cols, "metric")))) %>%
    dplyr::summarise(
      mean = round(mean(value, na.rm = TRUE), 3),
      sd = round(stats::sd(value, na.rm = TRUE), 3),
      n = sum(!is.na(value)),
      .groups = "drop"
    ) %>%
    dplyr::arrange(metric)
}

