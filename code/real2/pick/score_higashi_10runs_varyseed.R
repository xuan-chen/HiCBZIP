# score_higashi_10runs_varyseed.R
# Run with: Rscript score_higashi_10runs_varyseed.R

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

scool_path <- file.path("examples", "data", "oocyte_zygote_mm10_1M.scool")
ref_path <- file.path("examples", "data", "oocyte_zygote_ref")
if (!file.exists(scool_path)) stop("Missing scool: ", scool_path)
if (!file.exists(ref_path)) stop("Missing reference: ", ref_path)

n_runs <- 10L
seed_start <- 2026L
force_rerun <- TRUE
force_cpu <- TRUE
use_xy <- TRUE
no_viz <- TRUE

assembly <- "mm10"
resolution <- "1M"
min_depth <- "5000"
higashi_epochs <- 60L

embedding_alg <- "higashi"
embedding_dir <- "higashi"
embedding_json <- "higashi.json"
out_root <- file.path("results", "oocyte_zygote_higashi_10runs")
log_root <- file.path(out_root, "run_logs")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
dir.create(log_root, recursive = TRUE, showWarnings = FALSE)

metric_keys <- c(
  "ari_k-means", "ari_agglomerative", "ari_gmm", "ari_louvain", "ari_leiden", "best_ari",
  "nmi_k-means", "nmi_agglomerative", "nmi_gmm", "nmi_louvain", "nmi_leiden", "best_nmi",
  "silhouette_k-means", "silhouette_agglomerative", "silhouette_gmm", "silhouette_louvain",
  "silhouette_leiden", "best_silhouette", "best_silhouette-gt", "wall_time"
)

run_score_cmd <- function(score_bin, args, label, stdout_file, stderr_file, force_cpu = FALSE) {
  cat("\n[", label, "]\n", sep = "")
  if (force_cpu) cat("CUDA_VISIBLE_DEVICES=\"\" ")
  cat("score ", paste(args, collapse = " "), "\n", sep = "")

  old_cuda <- Sys.getenv("CUDA_VISIBLE_DEVICES", unset = NA_character_)
  if (force_cpu) Sys.setenv(CUDA_VISIBLE_DEVICES = "")
  on.exit({
    if (force_cpu) {
      if (is.na(old_cuda)) {
        Sys.unsetenv("CUDA_VISIBLE_DEVICES")
      } else {
        Sys.setenv(CUDA_VISIBLE_DEVICES = old_cuda)
      }
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

all_results <- vector("list", n_runs)

for (i in seq_len(n_runs)) {
  seed_i <- seed_start + i - 1L
  dset_name <- sprintf("oocyte_zygote_higashi_r%02d", i)
  metrics_json <- file.path(out_root, dset_name, resolution, embedding_dir, embedding_json)
  stdout_file <- file.path(log_root, sprintf("run_%02d_seed_%d.stdout.log", i, seed_i))
  stderr_file <- file.path(log_root, sprintf("run_%02d_seed_%d.stderr.log", i, seed_i))

  if (force_rerun || !file.exists(metrics_json)) {
    args <- c(
      "embed",
      "--dset", dset_name,
      "--out", out_root,
      "--assembly", assembly,
      "--resolution", resolution,
      "--scool", scool_path,
      "--reference", ref_path,
      "--embedding_algs", embedding_alg,
      "--min_depth", min_depth,
      "--higashi_epochs", as.character(higashi_epochs),
      "--seed", as.character(seed_i)
    )
    if (use_xy) args <- c(args, "--use_xy")
    if (no_viz) args <- c(args, "--no_viz")
    run_score_cmd(
      score_bin,
      args,
      label = paste("score embed higashi run", i),
      stdout_file = stdout_file,
      stderr_file = stderr_file,
      force_cpu = force_cpu
    )
  } else {
    cat("\n[skip] Reusing existing:", metrics_json, "\n")
  }

  out <- tibble(method = "Higashi", run = i, seed = seed_i)
  for (k in metric_keys) out[[k]] <- extract_metric(metrics_json, k)
  all_results[[i]] <- out
}

results_df <- bind_rows(all_results)
print(results_df)

results_long <- results_df %>%
  pivot_longer(cols = all_of(metric_keys), names_to = "metric", values_to = "value")

summary_df <- results_long %>%
  group_by(method, metric) %>%
  summarise(
    mean = round(mean(value, na.rm = TRUE), 3),
    sd = round(sd(value, na.rm = TRUE), 3),
    n = sum(!is.na(value)),
    .groups = "drop"
  ) %>%
  arrange(metric)

print(summary_df)

write.csv(results_df, file.path(out_root, "per_run_metrics.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_root, "summary_mean_sd.csv"), row.names = FALSE)

cat("Saved output tables to:", out_root, "\n")
cat("Saved per-run stdout/stderr logs to:", log_root, "\n")
