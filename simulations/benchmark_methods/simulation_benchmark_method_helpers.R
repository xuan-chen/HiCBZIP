# Shared helpers for simulation benchmark-method workflows.

simulation_chromosomes <- function() {
  c(1, 4, 5, 8, 10, 11, 15, 16, 17, 22)
}

simulation_coverages <- function() {
  c("0.01", "0.03", "0.05", "0.1", "0.2", "0.3", "0.4", "0.5", "0.7", "1")
}

simulation_region_info <- function() {
  data.frame(
    chr = simulation_chromosomes(),
    start_bp = c(50, 140, 100, 135, 25, 85, 80, 45, 55, 10) * 1e6,
    end_bp = c(55, 145, 105, 140, 30, 90, 85, 50, 60, 15) * 1e6
  ) |>
    dplyr::mutate(
      chr_name = paste0("chr", .data$chr),
      region = paste0(.data$chr_name, ":", .data$start_bp / 1e6, "M-", .data$end_bp / 1e6, "M")
    )
}

simulation_list_input_files <- function(input_dir) {
  pattern <- "^Simulation_snm3Cseq_human_brain_astrocytes_50k_chr[0-9XY]+_.*_K3X10\\.RData$"
  list.files(input_dir, pattern = pattern, full.names = TRUE)
}

simulation_extract_chr <- function(path) {
  as.integer(sub(".*chr([0-9XY]+).*", "\\1", basename(path)))
}

simulation_infer_n_bins_offdiag <- function(n_pairs) {
  n_bins <- as.integer(round((1 + sqrt(1 + 8 * n_pairs)) / 2))
  if (n_bins * (n_bins - 1) / 2 != n_pairs) {
    stop("Expected off-diagonal long-vector length n(n-1)/2, found ", n_pairs, call. = FALSE)
  }
  n_bins
}

simulation_make_pair_info_offdiag <- function(chr_name, n_bins, start_bp, resolution) {
  ij <- which(lower.tri(matrix(0, n_bins, n_bins), diag = FALSE), arr.ind = TRUE)
  data.frame(
    pair_id = seq_len(nrow(ij)),
    chrom1 = chr_name,
    pos1 = as.integer(start_bp + (ij[, 1] - 1L) * resolution),
    chrom2 = chr_name,
    pos2 = as.integer(start_bp + (ij[, 2] - 1L) * resolution)
  )
}

simulation_write_higashi_data_txt <- function(y_long_mat, cell_names, out_file, pair_info) {
  df <- as.data.frame(y_long_mat)
  colnames(df) <- cell_names
  df$pair_id <- seq_len(nrow(df))

  long_df <- tidyr::pivot_longer(df, cols = tidyselect::all_of(cell_names), names_to = "cell_name", values_to = "count")
  long_df <- dplyr::left_join(long_df, pair_info, by = "pair_id")
  long_df <- dplyr::mutate(
    long_df,
    cell_id = match(.data$cell_name, cell_names) - 1L,
    count = as.integer(.data$count)
  )
  long_df <- dplyr::select(long_df, tidyselect::all_of(c("cell_name", "cell_id", "chrom1", "pos1", "chrom2", "pos2", "count")))
  long_df <- dplyr::filter(long_df, .data$count > 0)
  write.table(long_df, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(long_df)
}

simulation_write_label_txt <- function(cell_names, out_file, default_label = 0L) {
  labels <- data.frame(
    cell_name = cell_names,
    cell_id = seq_along(cell_names) - 1L,
    label = rep(default_label, length(cell_names))
  )
  write.table(labels, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(labels)
}

simulation_read_csv_matrix <- function(path) {
  as.matrix(read.csv(path, header = FALSE, check.names = FALSE)) |>
    apply(2, as.numeric) |>
    as.matrix()
}
