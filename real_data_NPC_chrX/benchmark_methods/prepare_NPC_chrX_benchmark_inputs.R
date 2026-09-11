# Prepare per-coverage NPC chrX inputs for benchmark methods.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "tidyr", "jsonlite", "tidyselect"))
source_hicbzip_core()
source(file.path(script_dir, "..", "NPC_chrX_workflow_helpers.R"))

input_file <- Sys.getenv(
  "HICBZIP_REAL1_INPUT_RDATA",
  unset = path_here("data", "processed", "NPC_chrX", "data_NPC250k_0h_X_full.RData")
)
out_root <- Sys.getenv(
  "HICBZIP_REAL1_BENCHMARK_INPUT_DIR",
  unset = path_here("results", "NPC_chrX", "benchmark_inputs")
)

require_files(input_file, label = "NPC chrX processed input")
load(input_file)
if (!exists("true_muS")) true_muS <- replicate(K, bulk)

coverage_labels <- npc_chrX_coverage_labels()
coverage_values <- npc_chrX_coverage_values()
sim_y_list <- npc_chrX_make_sim_y_list(true_muS, coverage_labels = coverage_labels)
lambda_vals <- coverage_values
lambda_list <- paste0("ld=", coverage_labels)
resolution <- if (exists("resolution")) resolution else 250000L
n_bins <- if (exists("n")) n else npc_chrX_infer_n_bins_from_diag_long(nrow(true_muS))
k_cells <- ncol(true_muS)
cell_names <- colnames(true_muS)
if (is.null(cell_names)) cell_names <- paste0("cell", seq_len(k_cells))

ensure_dir(out_root)
save(sim_y_list, lambda_vals, lambda_list, coverage_labels, resolution, n_bins,
     file = file.path(out_root, "sim_y_list_NPC_chrX.RData"))

pair_info <- npc_chrX_make_pair_info_diag("chrX", n_bins = n_bins, start_bp = 0L, resolution = resolution)

sc_root <- file.path(out_root, "scHiCluster")
higashi_root <- file.path(out_root, "Higashi", "runs_chrX")
fast_root <- file.path(out_root, "FastHigashi", "runs_chrX")
for (path in c(sc_root, higashi_root, fast_root)) ensure_dir(path)

write_config <- function(run_dir, config_name, fast_higashi = FALSE) {
  cfg <- list(
    config_name = config_name,
    data_dir = normalizePath(run_dir, winslash = "/", mustWork = FALSE),
    temp_dir = normalizePath(file.path(run_dir, "temp"), winslash = "/", mustWork = FALSE),
    genome_reference_path = Sys.getenv("HICBZIP_HIGASHI_GENOME_REFERENCE", unset = "hg19.chrom.sizes.txt"),
    cytoband_path = Sys.getenv("HICBZIP_HIGASHI_CYTOBAND", unset = "cytoBand_hg19.txt"),
    chrom_list = list("chrX"),
    impute_list = list("chrX"),
    resolution = resolution,
    resolution_cell = resolution,
    minimum_distance = 0,
    maximum_distance = -1,
    minimum_impute_distance = 0,
    maximum_impute_distance = -1,
    local_transfer_range = 1,
    dimensions = 64,
    neighbor_num = 5,
    cpu_num = as.integer(Sys.getenv("HICBZIP_EXTERNAL_CPU", unset = "40")),
    gpu_num = as.integer(Sys.getenv("HICBZIP_EXTERNAL_GPU", unset = "1")),
    optional_smooth = FALSE,
    optional_quantile = FALSE,
    rank_thres = 1,
    loss_mode = "zinb",
    random_walk = FALSE,
    embedding_name = config_name
  )
  if (fast_higashi) {
    cfg$input_format <- "higashi_v1"
    cfg$structured <- TRUE
    cfg$resolution_fh <- list(resolution)
    cfg$batch_id <- "label"
  }
  writeLines(jsonlite::toJSON(cfg, pretty = TRUE, auto_unbox = TRUE), file.path(run_dir, "config.json"))
  invisible(cfg)
}

manifest_rows <- vector("list", length(sim_y_list))

for (ii in seq_along(sim_y_list)) {
  cov_label <- names(sim_y_list)[[ii]]
  cov_scientific <- npc_chrX_scientific_labels()[[ii]]
  y_i <- sim_y_list[[ii]]
  colnames(y_i) <- cell_names

  sc_dir <- file.path(sc_root, paste0("NPC250k_0h_X_full_", cov_label))
  sc_pair_dir <- file.path(sc_dir, "contact_pair")
  ensure_dir(sc_pair_dir)
  for (cell_idx in seq_len(ncol(y_i))) {
    nonzero <- which(y_i[, cell_idx] > 0)
    contact <- pair_info[nonzero, c("chrom1", "pos1", "chrom2", "pos2")]
    out_cell <- file.path(sc_pair_dir, paste0("NPC250k_0h_X_full_ds_", cov_label, "_", cell_idx, "_chrX.txt"))
    write.table(contact, out_cell, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }

  for (method_root in list(higashi_root, fast_root)) {
    is_fast <- identical(method_root, fast_root)
    run_dir <- file.path(method_root, paste0("chrX_cov_", cov_scientific))
    ensure_dir(run_dir)
    ensure_dir(file.path(run_dir, "temp"))
    ensure_dir(file.path(run_dir, "export_long_diag"))
    npc_chrX_write_higashi_data_txt(y_i, cell_names, file.path(run_dir, "data.txt"), pair_info)
    npc_chrX_write_label_txt(cell_names, file.path(run_dir, "label.txt"))
    write_config(run_dir, paste0(if (is_fast) "fasthigashi_" else "higashi_", "chrX_cov_", cov_scientific), fast_higashi = is_fast)
    metadata <- list(
      chr = "X",
      chr_name = "chrX",
      coverage = cov_scientific,
      repository_coverage = cov_label,
      resolution = resolution,
      start_bp = 0L,
      end_bp = as.integer(n_bins * resolution),
      n_bins = n_bins,
      vec_len = as.integer(n_bins * (n_bins + 1) / 2),
      n_cells = ncol(y_i),
      diagonal_included = TRUE
    )
    if (is_fast) {
      metadata$default_rank <- 3L
      metadata$conv_threshold <- 0.1
    }
    writeLines(jsonlite::toJSON(metadata, pretty = TRUE, auto_unbox = TRUE), file.path(run_dir, "run_metadata.json"))
  }

  manifest_rows[[ii]] <- data.frame(
    coverage = cov_label,
    scientific_coverage = cov_scientific,
    n_bins = n_bins,
    n_cells = ncol(y_i),
    vec_len = nrow(y_i)
  )
}

manifest <- dplyr::bind_rows(manifest_rows)
write.table(manifest, file.path(out_root, "NPC_chrX_benchmark_input_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved benchmark-method inputs under: ", normalizePath(out_root, winslash = "/", mustWork = FALSE))
