# Prepare simulation inputs for benchmark methods.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
require_packages(c("dplyr", "tidyr", "jsonlite", "tidyselect"))
source_hicbzip_core()
source(file.path(script_dir, "simulation_benchmark_method_helpers.R"))

input_dir <- Sys.getenv(
  "HICBZIP_SIM_INPUT_DIR",
  unset = path_here("data", "processed", "simulation", "input")
)
out_root <- Sys.getenv(
  "HICBZIP_SIM_EXTERNAL_INPUT_DIR",
  unset = path_here("results", "simulation", "benchmark_methods", "inputs")
)
resolution <- 50000L
coverage_labels <- simulation_coverages()
region_info <- simulation_region_info()

sim_files <- simulation_list_input_files(input_dir)
if (length(sim_files) == 0) stop("No processed simulation input files found in: ", input_dir, call. = FALSE)

for (path in c(out_root, file.path(out_root, "scHiCluster"), file.path(out_root, "Higashi", "runs"), file.path(out_root, "FastHigashi", "runs"))) {
  ensure_dir(path)
}

manifest_rows <- list()

for (input_file in sim_files) {
  env <- new.env(parent = emptyenv())
  load(input_file, envir = env)
  if (!exists("true_muS", envir = env)) stop("Missing true_muS in: ", input_file, call. = FALSE)
  true_muS <- as.matrix(env$true_muS)
  chr <- simulation_extract_chr(input_file)
  reg <- region_info[region_info$chr == chr, , drop = FALSE]
  if (nrow(reg) != 1) stop("No region metadata for chr", chr, call. = FALSE)

  n_bins <- simulation_infer_n_bins_offdiag(nrow(true_muS))
  pair_info <- simulation_make_pair_info_offdiag(reg$chr_name, n_bins, reg$start_bp, resolution)
  cell_names <- colnames(true_muS)
  if (is.null(cell_names)) cell_names <- paste0("cell", seq_len(ncol(true_muS)))

  for (coverage_label in coverage_labels) {
    set.seed(123456)
    sim_y <- matrix(
      rpois(nrow(true_muS) * ncol(true_muS), lambda = as.numeric(coverage_label) * as.numeric(true_muS)),
      nrow = nrow(true_muS),
      ncol = ncol(true_muS),
      dimnames = dimnames(true_muS)
    )

    sc_dir <- file.path(out_root, "scHiCluster", paste0("sim_HBA3_random_chr_", coverage_label))
    sc_pair_dir <- file.path(sc_dir, "contact_pair")
    ensure_dir(sc_pair_dir)
    for (cell_idx in seq_len(ncol(sim_y))) {
      nonzero <- which(sim_y[, cell_idx] > 0)
      contact <- pair_info[nonzero, c("chrom1", "pos1", "chrom2", "pos2")]
      out_cell <- file.path(sc_pair_dir, paste0("sim_HBA3_random_chr_ds_", coverage_label, "_", cell_idx, "_", reg$chr_name, ".txt"))
      write.table(contact, out_cell, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    }

    for (method_name in c("Higashi", "FastHigashi")) {
      run_dir <- file.path(out_root, method_name, "runs", paste0("chr", chr, "_cov", coverage_label))
      ensure_dir(run_dir)
      ensure_dir(file.path(run_dir, "temp"))
      ensure_dir(file.path(run_dir, "export_long"))
      simulation_write_higashi_data_txt(sim_y, cell_names, file.path(run_dir, "data.txt"), pair_info)
      simulation_write_label_txt(cell_names, file.path(run_dir, "label.txt"))

      cfg <- list(
        config_name = paste0(tolower(method_name), "_chr", chr, "_cov", coverage_label),
        data_dir = normalizePath(run_dir, winslash = "/", mustWork = FALSE),
        temp_dir = normalizePath(file.path(run_dir, "temp"), winslash = "/", mustWork = FALSE),
        genome_reference_path = Sys.getenv("HICBZIP_HIGASHI_GENOME_REFERENCE", unset = "hg19.chrom.sizes.txt"),
        cytoband_path = Sys.getenv("HICBZIP_HIGASHI_CYTOBAND", unset = "cytoBand_hg19.txt"),
        chrom_list = list(reg$chr_name),
        impute_list = list(reg$chr_name),
        resolution = resolution,
        resolution_cell = resolution,
        minimum_distance = 0,
        maximum_distance = -1,
        minimum_impute_distance = 0,
        maximum_impute_distance = -1,
        local_transfer_range = 1,
        dimensions = 64,
        neighbor_num = 5,
        cpu_num = as.integer(Sys.getenv("HICBZIP_EXTERNAL_CPU", unset = "8")),
        gpu_num = as.integer(Sys.getenv("HICBZIP_EXTERNAL_GPU", unset = "1")),
        optional_smooth = FALSE,
        optional_quantile = FALSE,
        rank_thres = 1,
        loss_mode = "zinb",
        random_walk = FALSE,
        embedding_name = paste0(tolower(method_name), "_chr", chr, "_cov", coverage_label)
      )
      if (method_name == "FastHigashi") {
        cfg$input_format <- "higashi_v1"
        cfg$structured <- TRUE
        cfg$resolution_fh <- list(resolution)
        cfg$batch_id <- "label"
      }
      writeLines(jsonlite::toJSON(cfg, pretty = TRUE, auto_unbox = TRUE), file.path(run_dir, "config.json"))
      metadata <- list(
        chr = chr,
        chr_name = reg$chr_name,
        coverage = coverage_label,
        region = reg$region,
        start_bp = as.integer(reg$start_bp),
        end_bp = as.integer(reg$end_bp),
        resolution = resolution,
        n_bins = n_bins,
        vec_len = as.integer(nrow(sim_y)),
        n_cells = as.integer(ncol(sim_y)),
        diagonal_included = FALSE,
        default_rank = 3L,
        conv_threshold = 0.1
      )
      writeLines(jsonlite::toJSON(metadata, pretty = TRUE, auto_unbox = TRUE), file.path(run_dir, "run_metadata.json"))
    }

    manifest_rows[[length(manifest_rows) + 1L]] <- data.frame(
      chr = chr,
      coverage = coverage_label,
      n_bins = n_bins,
      vec_len = nrow(sim_y),
      n_cells = ncol(sim_y)
    )
  }
}

manifest <- dplyr::bind_rows(manifest_rows) |>
  dplyr::arrange(.data$chr, as.numeric(.data$coverage))
write.table(manifest, file.path(out_root, "simulation_external_input_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
message("Saved simulation benchmark-method inputs under: ", normalizePath(out_root, winslash = "/", mustWork = FALSE))
