# Run SCORE Fast-Higashi for the oocyte-to-zygote benchmark.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()

source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
source(file.path(script_dir, "score_benchmark_helpers.R"))

score_bin <- resolve_score_bin()

raw_scool <- path_here("data", "bhzip_score_compare_from_pairs", "oocyte_zygote_raw_1M.scool")
ref_file <- path_here("data", "bhzip_score_compare_from_pairs", "oocyte_zygote_ref_min_depth_5000.tsv")
require_files(c(raw_scool, ref_file), label = "SCORE Fast-Higashi input")

n_runs <- 10L
seed_start <- 2026L
force_rerun <- TRUE
use_xy <- TRUE
no_viz <- TRUE

assembly <- "mm10"
resolution <- "1M"
min_depth <- "5000"

embedding_alg <- "fast_higashi"
embedding_dir <- "fast_higashi"
embedding_json <- "fast_higashi.json"
out_root <- path_here("results", "SCORE_oocyte_zygote", "oocyte_zygote_fast_higashi_10runs")
log_root <- file.path(out_root, "run_logs")
ensure_dir(out_root)
ensure_dir(log_root)

all_results <- vector("list", n_runs)

for (i in seq_len(n_runs)) {
  seed_i <- seed_start + i - 1L
  dset_name <- sprintf("oocyte_zygote_fast_higashi_r%02d", i)
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
      "--scool", raw_scool,
      "--reference", ref_file,
      "--embedding_algs", embedding_alg,
      "--min_depth", min_depth,
      "--seed", as.character(seed_i)
    )
    if (use_xy) args <- c(args, "--use_xy")
    if (no_viz) args <- c(args, "--no_viz")
    run_score_cmd(
      score_bin,
      args,
      label = paste("score embed Fast-Higashi run", i),
      stdout_file = stdout_file,
      stderr_file = stderr_file
    )
  } else {
    cat("\n[skip] Reusing existing:", metrics_json, "\n")
  }

  out <- tibble(method = "Fast-Higashi", run = i, seed = seed_i)
  for (k in score_metric_keys) out[[k]] <- extract_score_metric(metrics_json, k)
  all_results[[i]] <- out
}

results_df <- bind_rows(all_results)
summary_df <- summarize_score_results(results_df)

write.csv(results_df, file.path(out_root, "per_run_metrics.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_root, "summary_mean_sd.csv"), row.names = FALSE)

cat("Saved Fast-Higashi SCORE outputs to:", out_root, "\n")
cat("Saved per-run stdout/stderr logs to:", log_root, "\n")

