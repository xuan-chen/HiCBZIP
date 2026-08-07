# score_innerproduct_raw_vs_bhzip_nm_10runs_varyseed_plusviz.R
# Run with: Rscript score_innerproduct_raw_vs_bhzip_nm_10runs_varyseed_plusviz.R

required_pkgs <- c("jsonlite", "dplyr", "tibble", "tidyr", "ggplot2", "png")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))

library(jsonlite)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(png)
library(grid)

score_bin <- Sys.which("score")
if (!nzchar(score_bin)) {
  fallback <- c("/home/xchen/miniconda3/envs/score/bin/score", "C:/Users/67402/.conda/envs/score_py39/Scripts/score.exe")
  hit <- fallback[file.exists(fallback)]
  if (length(hit) > 0) score_bin <- hit[[1]]
}
if (!nzchar(score_bin)) stop("`score` command not found.")

root_dir <- file.path("data", "bhzip_score_compare_from_pairs_nm")
raw_scool <- file.path(root_dir, "oocyte_zygote_raw_1M.scool")
imp_scool <- file.path(root_dir, "oocyte_zygote_bhzip_nm_1M.scool")
ref_file <- file.path(root_dir, "oocyte_zygote_ref_min_depth_5000.tsv")

if (!file.exists(raw_scool)) stop("Missing raw scool: ", raw_scool)
if (!file.exists(imp_scool)) stop("Missing N(M) scool: ", imp_scool)
if (!file.exists(ref_file)) stop("Missing reference: ", ref_file)

n_runs <- 10L
seed_start <- 2026L
force_rerun <- TRUE
embedding_alg <- "InnerProduct"
embedding_dir <- "innerproduct"
embedding_json <- "innerproduct.json"

metric_keys <- c(
  "ari_k-means", "ari_agglomerative", "ari_gmm", "ari_louvain", "ari_leiden", "best_ari",
  "silhouette_k-means", "silhouette_agglomerative", "silhouette_gmm", "silhouette_louvain",
  "silhouette_leiden", "best_silhouette", "best_silhouette-gt"
)

run_score_cmd <- function(score_bin, args, label) {
  cat("\n[", label, "]\n", sep = "")
  cat("score ", paste(args, collapse = " "), "\n", sep = "")
  stdout_file <- tempfile(fileext = ".log")
  stderr_file <- tempfile(fileext = ".log")
  status <- system2(score_bin, args = args, stdout = stdout_file, stderr = stderr_file)
  out <- if (file.exists(stdout_file)) readLines(stdout_file, warn = FALSE) else character(0)
  err <- if (file.exists(stderr_file)) readLines(stderr_file, warn = FALSE) else character(0)
  cat("Return code:", status, "\n")
  if (length(out) > 0) cat(paste(out, collapse = "\n"), "\n")
  if (!identical(status, 0L)) {
    cat("\n--- STDERR ---\n")
    if (length(err) > 0) cat(paste(err, collapse = "\n"), "\n")
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

run_one <- function(method_name, scool_path, run_idx, seed) {
  dset_name <- sprintf("oocyte_zygote_%s_pairs_r%02d", method_name, run_idx)
  metrics_json <- file.path("results", dset_name, "1M", embedding_dir, embedding_json)

  if (force_rerun || !file.exists(metrics_json)) {
    args <- c(
      "embed",
      "--dset", dset_name,
      "--resolution", "1M",
      "--scool", scool_path,
      "--reference", ref_file,
      "--embedding_algs", embedding_alg,
      "--min_depth", "5000",
      "--use_xy",
      "--seed", as.character(seed)
    )
    run_score_cmd(score_bin, args, label = paste("score embed", method_name, "run", run_idx))
  }

  out <- tibble(method = method_name, run = run_idx, seed = seed)
  for (k in metric_keys) out[[k]] <- extract_metric(metrics_json, k)
  out
}

all_results <- list()
idx <- 1L
for (i in seq_len(n_runs)) {
  seed_i <- seed_start + i - 1L
  all_results[[idx]] <- run_one("raw_nm", raw_scool, i, seed_i); idx <- idx + 1L
  all_results[[idx]] <- run_one("bhzip_nm", imp_scool, i, seed_i); idx <- idx + 1L
}

results_df <- bind_rows(all_results)
results_long <- results_df %>% pivot_longer(cols = all_of(metric_keys), names_to = "metric", values_to = "value")
summary_df <- results_long %>% group_by(method, metric) %>% summarise(mean = round(mean(value, na.rm = TRUE),3), sd = round(sd(value, na.rm = TRUE),3), n = sum(!is.na(value)), .groups = "drop")

final_report_df <- summary_df %>%
  select(method, metric, mean, sd) %>%
  pivot_wider(names_from = method, values_from = c(mean, sd), names_glue = "{method}_{.value}") %>%
  mutate(
    raw_nm_mean_sd = sprintf("%.3f(%.3f)", raw_nm_mean, raw_nm_sd),
    bhzip_nm_mean_sd = sprintf("%.3f(%.3f)", bhzip_nm_mean, bhzip_nm_sd),
    delta = round(bhzip_nm_mean - raw_nm_mean, 3)
  ) %>%
  select(metric, raw_nm_mean_sd, bhzip_nm_mean_sd, delta) %>%
  arrange(desc(delta))

out_dir <- file.path("results", "oocyte_zygote_raw_vs_bhzip_nm_10runs_innerproduct")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results_df, file.path(out_dir, "per_run_metrics.csv"), row.names = FALSE)
write.csv(summary_df, file.path(out_dir, "summary_mean_sd.csv"), row.names = FALSE)
write.csv(final_report_df, file.path(out_dir, "final_report_mean_sd.csv"), row.names = FALSE)

# quick metric plots
plot_metrics <- c("best_ari", "best_silhouette", "best_silhouette-gt")
p_box <- ggplot(results_long %>% filter(metric %in% plot_metrics), aes(x = method, y = value, fill = method)) +
  geom_boxplot(width = 0.55, alpha = 0.65, outlier.shape = NA) +
  geom_jitter(width = 0.08, alpha = 0.8, size = 1.8) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c(raw_nm = "#4E79A7", bhzip_nm = "#E15759")) +
  labs(title = "InnerProduct Metrics: Raw vs HiCBZIP-N(M)", x = "Input", y = "Metric Value") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

ggsave(file.path(out_dir, "innerproduct_nm_raw_vs_bhzip_nm_boxplot_selected_metrics.png"), p_box, width = 10, height = 6, dpi = 200)

# run-1 side-by-side embedding panel
embedding_png_path <- function(method_name, run_idx) {
  dset_name <- sprintf("oocyte_zygote_%s_pairs_r%02d", method_name, run_idx)
  file.path("results", dset_name, "1M", embedding_dir, "celltype_plots", "embedding.png")
}
raw_png <- embedding_png_path("raw_nm", 1)
imp_png <- embedding_png_path("bhzip_nm", 1)
if (file.exists(raw_png) && file.exists(imp_png)) {
  raw_img <- readPNG(raw_png)
  imp_img <- readPNG(imp_png)
  png(file.path(out_dir, "run_01_seed_2026_raw_vs_bhzip_nm_embedding.png"), width = 2400, height = 1200, res = 180)
  grid.newpage()
  grid.text("InnerProduct Embedding Comparison - Run 01 (Seed 2026)", x = 0.5, y = 0.98, gp = gpar(fontsize = 18, fontface = "bold"))
  pushViewport(viewport(x = 0.5, y = 0.46, width = 0.98, height = 0.86, layout = grid.layout(nrow = 1, ncol = 2)))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1)); grid.raster(raw_img, interpolate = FALSE); grid.text("Raw", x=0.5, y=0.97, gp=gpar(fontsize=15,fontface="bold",col="#4E79A7")); popViewport()
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2)); grid.raster(imp_img, interpolate = FALSE); grid.text("HiCBZIP-N(M)", x=0.5, y=0.97, gp=gpar(fontsize=15,fontface="bold",col="#E15759")); popViewport()
  popViewport(); dev.off()
}

cat("Saved outputs to:", out_dir, "\n")
