# Run SCORE scHiCluster for the oocyte-to-zygote benchmark.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()

source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
source(file.path(script_dir, "score_benchmark_helpers.R"))

score_bin <- resolve_score_bin()

raw_scool <- path_here("data", "processed", "SCORE_oocyte_zygote", "score_ready_inputs", "HiCBZIP_GB_NB", "oocyte_zygote_raw_1M.scool")
ref_file <- path_here("data", "processed", "SCORE_oocyte_zygote", "score_ready_inputs", "HiCBZIP_GB_NB", "oocyte_zygote_ref_min_depth_5000.tsv")
require_files(c(raw_scool, ref_file), label = "SCORE scHiCluster input")

n_runs <- 10L
seed_start <- 2026L
force_rerun <- TRUE
use_xy <- TRUE
no_viz <- FALSE

embedding_alg <- "scHiCluster"
embedding_dir <- "schicluster:vc_sqrt_norm,convolution,random_walk"
embedding_json <- paste0(embedding_dir, ".json")
out_root <- path_here("results", "SCORE_oocyte_zygote", "oocyte_zygote_schicluster_10runs")
ensure_dir(out_root)

all_results <- vector("list", n_runs)

for (i in seq_len(n_runs)) {
  seed_i <- seed_start + i - 1L
  dset_name <- sprintf("oocyte_zygote_schicluster_r%02d", i)
  metrics_json <- file.path(out_root, dset_name, "1M", embedding_dir, embedding_json)

  if (force_rerun || !file.exists(metrics_json)) {
    args <- c(
      "embed",
      "--dset", dset_name,
      "--out", out_root,
      "--resolution", "1M",
      "--scool", raw_scool,
      "--reference", ref_file,
      "--embedding_algs", embedding_alg,
      "--min_depth", "5000",
      "--seed", as.character(seed_i)
    )
    if (use_xy) args <- c(args, "--use_xy")
    if (no_viz) args <- c(args, "--no_viz")
    run_score_cmd(score_bin, args, label = paste("score embed scHiCluster run", i))
  } else {
    cat("\n[skip] Reusing existing:", metrics_json, "\n")
  }

  out <- tibble(method = "scHiCluster", run = i, seed = seed_i)
  for (k in score_metric_keys) out[[k]] <- extract_score_metric(metrics_json, k)
  all_results[[i]] <- out
}

results_df <- bind_rows(all_results)
summary_df <- summarize_score_results(results_df)

write.csv(results_df, file.path(out_root, "per_run_metrics.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_root, "summary_mean_sd.csv"), row.names = FALSE)

cat("Saved scHiCluster SCORE outputs to:", out_root, "\n")
