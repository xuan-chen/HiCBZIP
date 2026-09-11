# Shared helpers for the NPC chrX manuscript workflow.

npc_chrX_coverage_labels <- function() {
  c("1e-04", "2e-04", "5e-04", "0.001", "0.002", "0.005", "0.01")
}

npc_chrX_coverage_values <- function() {
  as.numeric(npc_chrX_coverage_labels())
}

npc_chrX_scientific_labels <- function() {
  c("1e-04", "2e-04", "5e-04", "1e-03", "2e-03", "5e-03", "1e-02")
}

npc_chrX_normalize_cov_label <- function(x) {
  values <- suppressWarnings(as.numeric(as.character(x)))
  out <- rep(NA_character_, length(values))
  labels <- npc_chrX_coverage_labels()
  for (ii in seq_along(values)) {
    hit <- which(abs(values[[ii]] - as.numeric(labels)) < 1e-12)
    if (length(hit) == 1) out[[ii]] <- labels[[hit]]
  }
  if (anyNA(out)) {
    bad <- paste(as.character(x)[is.na(out)], collapse = ", ")
    stop("Unknown NPC chrX coverage label(s): ", bad, call. = FALSE)
  }
  out
}

npc_chrX_infer_n_bins_from_diag_long <- function(n_long) {
  n_bins <- (-1 + sqrt(1 + 8 * n_long)) / 2
  n_bins <- as.integer(round(n_bins))
  if (n_bins * (n_bins + 1) / 2 != n_long) {
    stop("Cannot infer number of bins from diagonal-included long vector length: ", n_long, call. = FALSE)
  }
  n_bins
}

npc_chrX_make_sim_y_list <- function(true_muS, coverage_labels = npc_chrX_coverage_labels(), seed = 123456) {
  out <- lapply(coverage_labels, function(coverage) {
    set.seed(seed)
    matrix(
      rpois(nrow(true_muS) * ncol(true_muS), lambda = as.numeric(coverage) * true_muS),
      nrow = nrow(true_muS),
      ncol = ncol(true_muS),
      dimnames = dimnames(true_muS)
    )
  })
  stats::setNames(out, coverage_labels)
}

npc_chrX_make_pair_info_diag <- function(chr_name, n_bins, start_bp, resolution) {
  ij <- which(lower.tri(matrix(0, n_bins, n_bins), diag = TRUE), arr.ind = TRUE)
  data.frame(
    pair_id = seq_len(nrow(ij)),
    chrom1 = chr_name,
    pos1 = as.integer(start_bp + (ij[, 1] - 1L) * resolution),
    chrom2 = chr_name,
    pos2 = as.integer(start_bp + (ij[, 2] - 1L) * resolution)
  )
}

npc_chrX_write_higashi_data_txt <- function(y_long_mat, cell_names, out_file, pair_info) {
  stopifnot(is.matrix(y_long_mat), nrow(y_long_mat) == nrow(pair_info), ncol(y_long_mat) == length(cell_names))
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

npc_chrX_write_label_txt <- function(cell_names, out_file, default_label = 0L) {
  labels <- data.frame(
    cell_name = cell_names,
    cell_id = seq_along(cell_names) - 1L,
    label = rep(default_label, length(cell_names))
  )
  write.table(labels, file = out_file, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(labels)
}

npc_chrX_read_csv_matrix <- function(path) {
  as.matrix(read.csv(path, header = FALSE, check.names = FALSE)) |>
    apply(2, as.numeric) |>
    as.matrix()
}

npc_chrX_standardize_method_list <- function(x, coverage_labels = npc_chrX_coverage_labels(), method_name = NULL) {
  if (!is.list(x)) stop("Expected a list of per-coverage matrices.", call. = FALSE)
  if (length(x) != length(coverage_labels)) {
    stop("Expected ", length(coverage_labels), " coverage matrices, found ", length(x), call. = FALSE)
  }
  names(x) <- if (is.null(method_name)) coverage_labels else paste0(method_name, ",ld=", coverage_labels)
  x
}

npc_chrX_build_scc_table <- function(list_muS_combine, true_muS, k_cells, lambda_list, resolution) {
  scc_list <- lapply(list_muS_combine, function(muS) {
    sapply(seq_len(k_cells), function(k) evaluate_mu(true_muS[, k], muS[, k], smooth = 0, res = resolution))[5, ]
  })
  do.call(rbind, scc_list) |>
    as.data.frame() |>
    tibble::rownames_to_column("method_label") |>
    tidyr::pivot_longer(-.data$method_label, names_to = "cell_id", values_to = "SCC") |>
    dplyr::mutate(
      lambda = paste0("ld=", sub("^.*ld=", "", .data$method_label)),
      method = sub(",ld=.*$", "", .data$method_label)
    )
}

npc_chrX_build_is_table <- function(list_muS_combine, true_muS, k_cells, lambda_list, resolution) {
  is_list <- lapply(list_muS_combine, function(muS) {
    get_IS_cor(true_muS, muS, n_cell = k_cells, win_size = 5, resolution = resolution)
  })
  do.call(rbind, is_list) |>
    as.data.frame() |>
    tibble::rownames_to_column("method_label") |>
    tidyr::pivot_longer(-.data$method_label, names_to = "cell_id", values_to = "IS") |>
    dplyr::mutate(
      lambda = paste0("ld=", sub("^.*ld=", "", .data$method_label)),
      method = sub(",ld=.*$", "", .data$method_label)
    )
}
