# score_innerproduct_four_inputs_10runs_varyseed.R
# Run with: Rscript score_innerproduct_four_inputs_10runs_varyseed.R

required_pkgs <- c("jsonlite", "dplyr", "tibble", "tidyr")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}

library(jsonlite)
library(dplyr)
library(tibble)
library(tidyr)

score_bin <- Sys.which("score")
if (!nzchar(score_bin)) {
  fallback_candidates <- c(
    "/home/xchen/miniconda3/envs/score/bin/score",
    "C:/Users/67402/.conda/envs/score_py39/Scripts/score.exe"
  )
  existing <- fallback_candidates[file.exists(fallback_candidates)]
  if (length(existing) > 0) score_bin <- existing[[1]]
}
if (!nzchar(score_bin)) stop("`score` command not found.")

inputs <- tibble::tribble(
  ~input_code,    ~input_label,       ~scool_path,                                                                 ~ref_path,
  "raw",          "Raw",              file.path("data", "bhzip_score_compare_from_pairs", "oocyte_zygote_raw_1M.scool"),          file.path("data", "bhzip_score_compare_from_pairs", "oocyte_zygote_ref_min_depth_5000.tsv"),
  "bhzip",        "HiCBZIP-NB(GS)",   file.path("data", "bhzip_score_compare_from_pairs", "oocyte_zygote_bhzip_1M.scool"),        file.path("data", "bhzip_score_compare_from_pairs", "oocyte_zygote_ref_min_depth_5000.tsv"),
  "bhzip_nm",     "HiCBZIP-N(M)",     file.path("data", "bhzip_score_compare_from_pairs_nm", "oocyte_zygote_bhzip_nm_1M.scool"),  file.path("data", "bhzip_score_compare_from_pairs_nm", "oocyte_zygote_ref_min_depth_5000.tsv"),
  "schicimpute",  "scHiCImpute",      file.path("data", "schicimpute_score_compare_from_pairs", "oocyte_zygote_schicimpute_1M.scool"), file.path("data", "schicimpute_score_compare_from_pairs", "oocyte_zygote_ref_min_depth_5000.tsv")
)

missing_inputs <- inputs %>%
  filter(!file.exists(scool_path) | !file.exists(ref_path))
if (nrow(missing_inputs) > 0) {
  stop(
    "Missing input files:\n",
    paste(
      apply(missing_inputs, 1, function(x) paste(x[["input_code"]], "=>", x[["scool_path"]], "|", x[["ref_path"]])),
      collapse = "\n"
    )
  )
}

n_runs <- 10L
seed_start <- 2026L
force_rerun <- TRUE
use_xy <- TRUE
no_viz <- FALSE

embedding_alg <- "InnerProduct"
embedding_dir <- "innerproduct"
embedding_json <- "innerproduct.json"
out_root <- file.path("results", "oocyte_zygote_four_inputs_10runs_innerproduct")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

metric_keys <- c(
  "ari_k-means", "ari_agglomerative", "ari_gmm", "ari_louvain", "ari_leiden", "best_ari",
  "nmi_k-means", "nmi_agglomerative", "nmi_gmm", "nmi_louvain", "nmi_leiden", "best_nmi",
  "silhouette_k-means", "silhouette_agglomerative", "silhouette_gmm", "silhouette_louvain",
  "silhouette_leiden", "best_silhouette", "best_silhouette-gt", "wall_time"
)

run_score_cmd <- function(score_bin, args, label) {
  cat("\n[", label, "]\n", sep = "")
  cat("score ", paste(args, collapse = " "), "\n", sep = "")

  stdout_file <- tempfile(fileext = ".log")
  stderr_file <- tempfile(fileext = ".log")
  status <- system2(score_bin, args = args, stdout = stdout_file, stderr = stderr_file)

  stdout_lines <- if (file.exists(stdout_file)) readLines(stdout_file, warn = FALSE) else character(0)
  stderr_lines <- if (file.exists(stderr_file)) readLines(stderr_file, warn = FALSE) else character(0)

  cat("Return code:", status, "\n")
  if (length(stdout_lines) > 0) cat(paste(stdout_lines, collapse = "\n"), "\n")
  if (!identical(status, 0L)) {
    cat("\n--- STDERR ---\n")
    if (length(stderr_lines) > 0) cat(paste(stderr_lines, collapse = "\n"), "\n")
    stop(label, " failed.")
  }
}

extract_metric <- function(path, key) {
  if (!file.exists(path)) return(NA_real_)
  x <- fromJSON(path, simplifyVector = FALSE)[[key]]
  if (is.null(x)) return(NA_real_)
  if (is.list(x)) x <- unlist(x, recursive = TRUE, use.names = FALSE)
  v <- suppressWarnings(as.numeric(x))
  if (any(!is.na(v))) return(v[which(!is.na(v))[1]])
  NA_real_
}

run_one_experiment <- function(input_code, input_label, scool_path, ref_path, run_idx, seed) {
  dset_name <- sprintf("oocyte_zygote_%s_pairs_r%02d", input_code, run_idx)
  metrics_json <- file.path(out_root, dset_name, "1M", embedding_dir, embedding_json)

  if (force_rerun || !file.exists(metrics_json)) {
    args <- c(
      "embed",
      "--dset", dset_name,
      "--out", out_root,
      "--resolution", "1M",
      "--scool", scool_path,
      "--reference", ref_path,
      "--embedding_algs", embedding_alg,
      "--min_depth", "5000",
      "--seed", as.character(seed)
    )
    if (use_xy) args <- c(args, "--use_xy")
    if (no_viz) args <- c(args, "--no_viz")
    run_score_cmd(score_bin, args, label = paste("score embed", input_code, "run", run_idx))
  } else {
    cat("\n[skip] Reusing existing:", metrics_json, "\n")
  }

  out <- tibble(input_code = input_code, input_label = input_label, run = run_idx, seed = seed)
  for (k in metric_keys) out[[k]] <- extract_metric(metrics_json, k)
  out
}

all_results <- list()
idx <- 1L

for (i in seq_len(n_runs)) {
  seed_i <- seed_start + i - 1L
  for (j in seq_len(nrow(inputs))) {
    all_results[[idx]] <- run_one_experiment(
      input_code = inputs$input_code[[j]],
      input_label = inputs$input_label[[j]],
      scool_path = inputs$scool_path[[j]],
      ref_path = inputs$ref_path[[j]],
      run_idx = i,
      seed = seed_i
    )
    idx <- idx + 1L
  }
}

results_df <- bind_rows(all_results)
print(results_df)

results_long <- results_df %>%
  pivot_longer(cols = all_of(metric_keys), names_to = "metric", values_to = "value")

summary_df <- results_long %>%
  group_by(input_code, input_label, metric) %>%
  summarise(
    mean = round(mean(value, na.rm = TRUE), 3),
    sd = round(sd(value, na.rm = TRUE), 3),
    n = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  arrange(metric, input_code)

print(summary_df)

final_report_df <- summary_df %>%
  select(input_code, metric, mean, sd) %>%
  pivot_wider(
    names_from = input_code,
    values_from = c(mean, sd),
    names_glue = "{input_code}_{.value}"
  ) %>%
  mutate(
    raw_mean_sd = sprintf("%.3f(%.3f)", raw_mean, raw_sd),
    bhzip_mean_sd = sprintf("%.3f(%.3f)", bhzip_mean, bhzip_sd),
    bhzip_nm_mean_sd = sprintf("%.3f(%.3f)", bhzip_nm_mean, bhzip_nm_sd),
    schicimpute_mean_sd = sprintf("%.3f(%.3f)", schicimpute_mean, schicimpute_sd),
    gb_minus_raw = round(bhzip_mean - raw_mean, 3),
    nm_minus_raw = round(bhzip_nm_mean - raw_mean, 3),
    schicimpute_minus_raw = round(schicimpute_mean - raw_mean, 3)
  ) %>%
  select(metric, raw_mean_sd, bhzip_mean_sd, bhzip_nm_mean_sd, schicimpute_mean_sd,
         gb_minus_raw, nm_minus_raw, schicimpute_minus_raw)

print(final_report_df)

write.csv(results_df, file.path(out_root, "per_run_metrics.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_root, "summary_mean_sd.csv"), row.names = FALSE)
write.csv(final_report_df, file.path(out_root, "final_report_mean_sd.csv"), row.names = FALSE)

cat("Saved output tables to:", out_root, "\n")
