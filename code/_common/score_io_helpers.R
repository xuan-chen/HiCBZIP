ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

download_gse305523 <- function(data_dir = "data", force = FALSE) {
  ensure_dir(data_dir)
  raw_dir <- file.path(data_dir, "raw")
  ensure_dir(raw_dir)

  tar_urls <- c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305523/suppl/GSE305523_reanalysis_GSE84920_ramani.tar.gz",
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE305523&format=file&file=GSE305523_reanalysis_GSE84920_ramani.tar.gz"
  )
  meta_urls <- c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE305nnn/GSE305523/suppl/GSE305523_reanalysis_GSE84920_ramani_meta.tsv.gz",
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE305523&format=file&file=GSE305523_reanalysis_GSE84920_ramani_meta.tsv.gz"
  )

  tar_path <- file.path(raw_dir, "GSE305523_reanalysis_GSE84920_ramani.tar.gz")
  meta_path <- file.path(raw_dir, "GSE305523_reanalysis_GSE84920_ramani_meta.tsv.gz")
  extract_dir <- file.path(raw_dir, "GSE305523_reanalysis_GSE84920_ramani")
  ensure_dir(extract_dir)

  download_with_fallback <- function(urls, destfile) {
    last_err <- NULL
    for (u in urls) {
      ok <- tryCatch({
        utils::download.file(u, destfile = destfile, mode = "wb", quiet = FALSE)
        TRUE
      }, error = function(e) {
        last_err <<- conditionMessage(e)
        FALSE
      })
      if (ok && file.exists(destfile) && file.info(destfile)$size > 0) return(invisible(TRUE))
    }
    stop("All download attempts failed for ", basename(destfile), ". Last error: ", last_err)
  }

  if (force || !file.exists(tar_path) || file.info(tar_path)$size <= 0) {
    download_with_fallback(tar_urls, tar_path)
  }
  if (force || !file.exists(meta_path) || file.info(meta_path)$size <= 0) {
    download_with_fallback(meta_urls, meta_path)
  }
  already_extracted <- length(list.files(extract_dir, all.files = TRUE, no.. = TRUE)) > 0
  if (force || !already_extracted) {
    utils::untar(tar_path, exdir = extract_dir)
  }

  list(
    tar_path = normalizePath(tar_path, winslash = "/", mustWork = FALSE),
    metadata_path = normalizePath(meta_path, winslash = "/", mustWork = FALSE),
    extract_dir = normalizePath(extract_dir, winslash = "/", mustWork = FALSE)
  )
}
read_metadata_table <- function(path) {
  if (!file.exists(path)) stop("Metadata file not found: ", path)
  readr::read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

choose_first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) return(NA_character_)
  hit[[1]]
}

infer_metadata_columns <- function(meta_df) {
  cell_candidates <- c(
    "cell", "cell_id", "cellid", "Cell", "CellID", "cell_name", "Cellname",
    "barcode", "cell_barcode"
  )
  label_candidates <- c(
    "cell_type", "CellType", "celltype", "Cell type", "label",
    "cluster", "group", "stage", "condition"
  )

  cell_col <- choose_first_existing_col(meta_df, cell_candidates)
  label_col <- choose_first_existing_col(meta_df, label_candidates)

  if (is.na(cell_col)) {
    stop(
      "Could not infer cell ID column from metadata. Available columns:\n",
      paste(names(meta_df), collapse = ", ")
    )
  }
  if (is.na(label_col)) {
    stop(
      "Could not infer label column from metadata. Available columns:\n",
      paste(names(meta_df), collapse = ", ")
    )
  }

  list(cell_col = cell_col, label_col = label_col)
}

list_candidate_matrix_files <- function(extract_dir) {
  if (!dir.exists(extract_dir)) stop("Extract directory not found: ", extract_dir)

  patt <- "\\.(rds|mtx|mtx\\.gz|csv|csv\\.gz|tsv|tsv\\.gz|scool|cool|h5|hdf5)$"
  files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
  files <- files[grepl(patt, files, ignore.case = TRUE)]

  is_matrix_like <- grepl(
    "matrix|contact|count|hic|ramani|feature|cool",
    basename(files),
    ignore.case = TRUE
  )
  files2 <- files[is_matrix_like]
  if (length(files2) == 0) files2 <- files
  sort(unique(files2))
}

pick_default_matrix_file <- function(candidates) {
  if (length(candidates) == 0) stop("No candidate matrix-like files found.")

  score_file <- function(x) {
    b <- tolower(basename(x))
    score <- 0L
    if (grepl("ramani", b)) score <- score + 4L
    if (grepl("raw", b)) score <- score + 3L
    if (grepl("matrix|count|contact", b)) score <- score + 2L
    if (grepl("\\.rds$", b)) score <- score + 3L
    if (grepl("\\.mtx(\\.gz)?$", b)) score <- score + 2L
    score
  }

  scores <- vapply(candidates, score_file, integer(1))
  candidates <- candidates[order(-scores, candidates)]
  candidates[[1]]
}

read_vector_file <- function(path) {
  if (!file.exists(path)) stop("Vector file not found: ", path)
  x <- readr::read_lines(path, progress = FALSE)
  x <- x[nzchar(x)]
  x
}

read_contact_matrix <- function(matrix_path, cell_ids_path = NULL, feature_ids_path = NULL) {
  if (!file.exists(matrix_path)) stop("Matrix file not found: ", matrix_path)
  lower <- tolower(matrix_path)

  if (grepl("\\.rds$", lower)) {
    obj <- readRDS(matrix_path)

    if (inherits(obj, "Matrix") || is.matrix(obj)) {
      mat <- obj
      cell_ids <- colnames(mat)
      feature_ids <- rownames(mat)
    } else if (is.list(obj)) {
      mat <- obj$counts
      if (is.null(mat)) mat <- obj$matrix
      if (is.null(mat)) mat <- obj$X
      if (is.null(mat)) stop("RDS list does not contain counts/matrix/X.")
      cell_ids <- obj$cell_ids
      if (is.null(cell_ids)) cell_ids <- obj$cells
      if (is.null(cell_ids)) cell_ids <- colnames(mat)
      feature_ids <- obj$feature_ids
      if (is.null(feature_ids)) feature_ids <- obj$features
      if (is.null(feature_ids)) feature_ids <- rownames(mat)
    } else {
      stop("Unsupported RDS object class: ", class(obj)[1])
    }
  } else if (grepl("\\.mtx(\\.gz)?$", lower)) {
    mat <- Matrix::readMM(matrix_path)
    mat <- as(mat, "dgCMatrix")
    cell_ids <- NULL
    feature_ids <- NULL
  } else if (grepl("\\.csv(\\.gz)?$", lower)) {
    df <- readr::read_csv(matrix_path, show_col_types = FALSE, progress = FALSE)
    if (!is.numeric(df[[1]])) {
      feature_ids <- df[[1]]
      df <- df[, -1, drop = FALSE]
    } else {
      feature_ids <- NULL
    }
    mat <- as.matrix(df)
    cell_ids <- colnames(df)
  } else if (grepl("\\.tsv(\\.gz)?$", lower)) {
    df <- readr::read_tsv(matrix_path, show_col_types = FALSE, progress = FALSE)
    if (!is.numeric(df[[1]])) {
      feature_ids <- df[[1]]
      df <- df[, -1, drop = FALSE]
    } else {
      feature_ids <- NULL
    }
    mat <- as.matrix(df)
    cell_ids <- colnames(df)
  } else {
    stop("Unsupported matrix extension for: ", matrix_path)
  }

  if (!is.null(cell_ids_path)) {
    cell_ids <- read_vector_file(cell_ids_path)
  }
  if (!is.null(feature_ids_path)) {
    feature_ids <- read_vector_file(feature_ids_path)
  }

  if (is.null(cell_ids)) cell_ids <- colnames(mat)
  if (is.null(feature_ids)) feature_ids <- rownames(mat)

  list(mat = mat, cell_ids = cell_ids, feature_ids = feature_ids)
}

find_sidecar_ids <- function(matrix_path) {
  dir0 <- dirname(matrix_path)
  all_files <- list.files(dir0, full.names = TRUE)
  b <- basename(all_files)

  pick <- function(patterns) {
    hit <- all_files[vapply(
      patterns,
      function(p) grepl(p, b, ignore.case = TRUE),
      logical(length(b))
    )]
    if (length(hit) == 0) return(NULL)
    sort(unique(hit))[[1]]
  }

  list(
    cell_ids_path = pick(c("barcodes\\.tsv(\\.gz)?$", "cells\\.txt$", "cell_ids\\.txt$")),
    feature_ids_path = pick(c("features\\.tsv(\\.gz)?$", "genes\\.tsv(\\.gz)?$", "pairs\\.txt$", "feature_ids\\.txt$"))
  )
}

as_feature_by_cell <- function(mat, cell_ids, metadata_cells) {
  if (!is.numeric(mat) && !inherits(mat, "Matrix")) {
    stop("Matrix must be numeric/sparse numeric.")
  }
  if (any(mat < 0, na.rm = TRUE)) stop("Matrix has negative counts.")

  n_meta <- length(metadata_cells)
  by_col <- !is.null(cell_ids) && (ncol(mat) == length(cell_ids))
  by_row <- !is.null(cell_ids) && (nrow(mat) == length(cell_ids))

  if (by_col && sum(cell_ids %in% metadata_cells) >= max(10, floor(0.1 * n_meta))) {
    cell_axis <- "col"
  } else if (by_row && sum(cell_ids %in% metadata_cells) >= max(10, floor(0.1 * n_meta))) {
    mat <- Matrix::t(mat)
    cell_axis <- "row"
  } else {
    stop("Could not determine cell axis from matrix and metadata overlap.")
  }

  cell_ids2 <- if (cell_axis == "col") cell_ids else cell_ids
  if (length(cell_ids2) != ncol(mat)) {
    stop("Cell ID length does not match matrix cell dimension.")
  }
  colnames(mat) <- cell_ids2
  mat
}

align_matrix_and_metadata <- function(feature_by_cell, meta_df, cell_col, label_col) {
  keep <- !is.na(meta_df[[cell_col]]) & !is.na(meta_df[[label_col]])
  meta_df <- meta_df[keep, , drop = FALSE]

  common <- intersect(colnames(feature_by_cell), meta_df[[cell_col]])
  if (length(common) < 50) {
    stop("Too few matched cells between matrix and metadata: ", length(common))
  }

  meta_aligned <- dplyr::filter(meta_df, .data[[cell_col]] %in% common)
  meta_aligned <- dplyr::arrange(meta_aligned, match(.data[[cell_col]], common))
  mat_aligned <- feature_by_cell[, common, drop = FALSE]

  list(
    mat = mat_aligned,
    meta = meta_aligned
  )
}

stratified_sample_meta <- function(meta_df, label_col, n_target = 300, seed = 12345) {
  set.seed(seed)
  meta_df <- dplyr::mutate(meta_df, .label_tmp = as.character(.data[[label_col]]))
  meta_df <- meta_df[!is.na(meta_df$.label_tmp) & nzchar(meta_df$.label_tmp), , drop = FALSE]

  if (n_target >= nrow(meta_df)) return(meta_df)

  groups <- split(meta_df, meta_df$.label_tmp, drop = TRUE)
  labels <- names(groups)
  counts <- vapply(groups, nrow, integer(1))

  prop <- counts / sum(counts)
  n_per <- pmax(1L, floor(prop * n_target))
  names(n_per) <- labels

  # If we overshoot due to pmax(1), trim from largest groups first.
  while (sum(n_per) > n_target) {
    can_drop <- labels[n_per > 1L]
    if (length(can_drop) == 0) break
    key <- can_drop[[which.max(n_per[can_drop])]]
    n_per[key] <- n_per[key] - 1L
  }

  # If still short, add to largest groups cyclically.
  leftover <- n_target - sum(n_per)
  if (leftover > 0) {
    order_lbl <- labels[order(prop[labels], decreasing = TRUE)]
    for (i in seq_len(leftover)) {
      key <- order_lbl[[((i - 1) %% length(order_lbl)) + 1]]
      n_per[key] <- n_per[key] + 1L
    }
  }

  sampled_list <- lapply(labels, function(lbl) {
    df <- groups[[lbl]]
    take <- n_per[[lbl]]
    if (is.null(take) || is.na(take)) take <- 0L
    n_take <- min(nrow(df), as.integer(take))
    if (n_take <= 0L || nrow(df) == 0L) return(df[0, , drop = FALSE])
    df[sample.int(nrow(df), n_take), , drop = FALSE]
  })

  sampled <- dplyr::bind_rows(sampled_list)
  sampled <- dplyr::arrange(sampled, .data[[label_col]])
  sampled$.label_tmp <- NULL
  sampled
}

compute_bulk_referenced_coverage <- function(Y) {
  C_k <- Matrix::colSums(Y)
  C_bulk <- sum(C_k)
  K <- length(C_k)
  C_bulk_per_cell <- C_bulk / K
  coverage <- as.numeric(C_k / C_bulk)

  list(
    C_k = C_k,
    C_bulk = C_bulk,
    C_bulk_per_cell = C_bulk_per_cell,
    coverage = coverage
  )
}







upper_tri_rowmajor_id <- function(i, j, n) {
  # i,j are 1-based with i <= j
  before <- (i - 1) * (2 * n - i + 2) / 2
  as.integer(before + (j - i + 1))
}

list_scool_cells <- function(scool_path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required for .scool input. Install with install.packages('hdf5r').")
  }
  f <- hdf5r::H5File$new(scool_path, mode = "r")
  on.exit(try(f$close_all(), silent = TRUE), add = TRUE)

  if (!"cells" %in% names(f)) stop(".scool file does not contain '/cells' group: ", scool_path)
  names(f[["cells"]])
}

load_scool_feature_by_cell <- function(scool_path, cell_ids, target_chromosomes = "chr19", include_diag = TRUE, verbose = TRUE) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required for .scool input. Install with install.packages('hdf5r').")
  }
  f <- hdf5r::H5File$new(scool_path, mode = "r")
  on.exit(try(f$close_all(), silent = TRUE), add = TRUE)

  if (!"cells" %in% names(f)) stop(".scool file missing '/cells' group")
  if (!"bins" %in% names(f)) stop(".scool file missing '/bins' group")
  if (!"chroms" %in% names(f)) stop(".scool file missing '/chroms' group")

  all_cells <- names(f[["cells"]])
  keep_cells <- intersect(cell_ids, all_cells)
  if (length(keep_cells) == 0) stop("No overlapping cell IDs between metadata and .scool cells")

  chrom_names <- f[["chroms"]][["name"]][]
  chrom_names <- as.character(chrom_names)
  bins_chrom_idx <- as.integer(f[["bins"]][["chrom"]][]) + 1L
  if (any(bins_chrom_idx < 1 | bins_chrom_idx > length(chrom_names))) {
    stop("Invalid chrom indices in .scool bins/chrom")
  }
  bins_chr <- chrom_names[bins_chrom_idx]

  if (length(target_chromosomes) == 1 && is.na(target_chromosomes)) {
    use_bins <- seq_along(bins_chr)
  } else {
    use_bins <- which(bins_chr %in% target_chromosomes)
  }
  if (length(use_bins) == 0) {
    stop("No bins found for target chromosomes: ", paste(target_chromosomes, collapse = ","))
  }

  # Keep a contiguous local index over selected bins.
  old2new <- integer(length(bins_chr))
  old2new[use_bins] <- seq_along(use_bins)
  n_bins <- length(use_bins)

  if (!include_diag) {
    stop("include_diag = FALSE is not implemented for .scool loader in this notebook")
  }
  n_pairs <- as.integer(n_bins * (n_bins + 1) / 2)
  Y <- matrix(0, nrow = n_pairs, ncol = length(keep_cells))
  colnames(Y) <- keep_cells

  for (k in seq_along(keep_cells)) {
    cid <- keep_cells[[k]]
    g <- f[["cells"]][[cid]]
    if (!"pixels" %in% names(g)) next

    b1 <- as.integer(g[["pixels"]][["bin1_id"]][]) + 1L
    b2 <- as.integer(g[["pixels"]][["bin2_id"]][]) + 1L
    cnt <- as.numeric(g[["pixels"]][["count"]][])

    in_chr <- old2new[b1] > 0 & old2new[b2] > 0
    if (!any(in_chr)) next

    i <- old2new[b1[in_chr]]
    j <- old2new[b2[in_chr]]
    cts <- cnt[in_chr]

    # Ensure upper-triangle indexing i <= j for run_BZIP_GB_NB neighbor map.
    ii <- pmin(i, j)
    jj <- pmax(i, j)
    ids <- upper_tri_rowmajor_id(ii, jj, n_bins)

    # Aggregate duplicate pair ids.
    s <- tapply(cts, ids, sum)
    Y[as.integer(names(s)), k] <- as.numeric(s)

    if (verbose && (k %% 50 == 0 || k == length(keep_cells))) {
      message("Loaded .scool cell ", k, "/", length(keep_cells))
    }
  }

  rownames(Y) <- paste0("pair_", seq_len(n_pairs))
  Y
}


normalize_cell_id <- function(x) {
  x <- as.character(x)
  x <- tolower(x)
  gsub("[^a-z0-9]", "", x)
}

extract_barcode_keys <- function(x) {
  x0 <- as.character(x)
  x0 <- gsub("\\.1m$", "", x0, ignore.case = TRUE)
  x0 <- gsub("\\.txt$", "", x0, ignore.case = TRUE)
  after_us <- sub("^.*_", "", x0)
  pair <- toupper(after_us)

  second <- ifelse(grepl("-", pair), sub("^.*-", "", pair), pair)
  data.frame(
    .pair = pair,
    .second = second,
    stringsAsFactors = FALSE
  )
}

match_metadata_to_scool_cells <- function(meta_df, cell_col, label_col, scool_cells, min_cells = 50) {
  meta0 <- meta_df %>%
    dplyr::filter(!is.na(.data[[cell_col]]), !is.na(.data[[label_col]]))

  # 1) Direct exact match.
  direct <- meta0 %>%
    dplyr::filter(.data[[cell_col]] %in% scool_cells) %>%
    dplyr::mutate(.scool_cell = .data[[cell_col]], .match_method = "direct")

  # 2) Normalized full-string match.
  norm_scool <- tibble::tibble(
    .scool_cell = scool_cells,
    .norm = normalize_cell_id(scool_cells)
  ) %>%
    dplyr::add_count(.norm, name = "n_scool_norm") %>%
    dplyr::filter(n_scool_norm == 1)

  norm_meta <- meta0 %>%
    dplyr::mutate(.norm = normalize_cell_id(.data[[cell_col]])) %>%
    dplyr::add_count(.norm, name = "n_meta_norm") %>%
    dplyr::filter(n_meta_norm == 1)

  norm_map <- norm_meta %>%
    dplyr::inner_join(norm_scool %>% dplyr::select(.scool_cell, .norm), by = ".norm") %>%
    dplyr::mutate(.match_method = "normalized")

  # 3) Barcode-pair match (token after underscore, strip .txt/.1M).
  scool_keys <- cbind(
    tibble::tibble(.scool_cell = scool_cells),
    extract_barcode_keys(scool_cells)
  )
  meta_keys <- cbind(
    meta0,
    extract_barcode_keys(meta0[[cell_col]])
  )

  pair_map <- meta_keys %>%
    dplyr::add_count(.pair, name = "n_meta_pair") %>%
    dplyr::filter(n_meta_pair == 1, nzchar(.pair)) %>%
    dplyr::inner_join(
      scool_keys %>% dplyr::add_count(.pair, name = "n_scool_pair") %>% dplyr::filter(n_scool_pair == 1) %>% dplyr::select(.scool_cell, .pair),
      by = ".pair"
    ) %>%
    dplyr::mutate(.match_method = "barcode_pair")

  # 4) Fallback: second barcode token after '-'.
  second_map <- meta_keys %>%
    dplyr::add_count(.second, name = "n_meta_second") %>%
    dplyr::filter(n_meta_second == 1, nzchar(.second)) %>%
    dplyr::inner_join(
      scool_keys %>% dplyr::add_count(.second, name = "n_scool_second") %>% dplyr::filter(n_scool_second == 1) %>% dplyr::select(.scool_cell, .second),
      by = ".second"
    ) %>%
    dplyr::mutate(.match_method = "barcode_second")

  mapped <- dplyr::bind_rows(direct, norm_map, pair_map, second_map) %>%
    dplyr::distinct(.scool_cell, .keep_all = TRUE) %>%
    dplyr::arrange(match(.scool_cell, scool_cells))

  if (nrow(mapped) < min_cells) {
    ex_meta <- paste(head(as.character(meta0[[cell_col]]), 5), collapse = ", ")
    ex_scool <- paste(head(as.character(scool_cells), 5), collapse = ", ")
    stop(
      "Too few matched cells between metadata and .scool after direct/normalized/barcode matching: ", nrow(mapped),
      "\nExample metadata IDs: ", ex_meta,
      "\nExample .scool IDs: ", ex_scool
    )
  }

  mapped[[cell_col]] <- mapped$.scool_cell
  mapped
}



list_scool_chromosomes <- function(scool_path) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package 'hdf5r' is required for .scool input. Install with install.packages('hdf5r').")
  }
  f <- hdf5r::H5File$new(scool_path, mode = "r")
  on.exit(try(f$close_all(), silent = TRUE), add = TRUE)
  if (!"chroms" %in% names(f)) stop(".scool file does not contain '/chroms' group: ", scool_path)
  as.character(f[["chroms"]][["name"]][])
}

