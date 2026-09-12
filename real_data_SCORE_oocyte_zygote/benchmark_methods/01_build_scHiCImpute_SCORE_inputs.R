# Build scHiCImpute SCORE-ready pair files and .scool input for real-data study 2.

script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()

source(file.path(script_dir, "..", "..", "_common", "project_paths.R"))
source(file.path(script_dir, "score_benchmark_helpers.R"))
source_hicbzip_core(include_score_helpers = TRUE)

required_pkgs <- c("HiCImpute", "Matrix")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

score_bin <- resolve_score_bin()

raw_data_dir <- path_here("data", "processed", "SCORE_oocyte_zygote", "oocyte_zygote_mm10", "1M")
anchor_file <- path_here("data", "processed", "SCORE_oocyte_zygote", "mm10.genome_split_1M")
ref_file <- path_here("data", "processed", "SCORE_oocyte_zygote", "oocyte_zygote_ref")

require_dirs(raw_data_dir, label = "SCORE raw pair directory")
require_files(c(anchor_file, ref_file), label = "SCORE processed input")

cfg <- list(
  min_depth = 5000L,
  mc.cores = as.integer(Sys.getenv("HICBZIP_HICIMPUTE_CORES", unset = "8")),
  cutoff = 0.5,
  niter = 5000L,
  burnin = 500L,
  chr_subset = character(0)
)

normalize_key <- function(x) {
  x <- tolower(as.character(x))
  gsub("\\.1m$", "", x, ignore.case = TRUE)
}

read_pair_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(bin1 = integer(0), bin2 = integer(0), count = numeric(0)))
  }
  x <- tryCatch(read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE), error = function(e) data.frame())
  if (nrow(x) == 0 || ncol(x) < 3) {
    return(data.frame(bin1 = integer(0), bin2 = integer(0), count = numeric(0)))
  }
  x <- x[, 1:3]
  colnames(x) <- c("bin1", "bin2", "count")
  x$bin1 <- as.integer(x$bin1)
  x$bin2 <- as.integer(x$bin2)
  x$count <- as.numeric(x$count)
  x[!is.na(x$bin1) & !is.na(x$bin2) & !is.na(x$count), , drop = FALSE]
}

build_chr_long_offdiag <- function(dt, chr_ids) {
  n_bins <- length(chr_ids)
  if (n_bins < 2) return(numeric(0))

  old2new <- integer(max(chr_ids) + 1L)
  old2new[chr_ids + 1L] <- seq_len(n_bins)
  in_chr <- dt$bin1 %in% chr_ids & dt$bin2 %in% chr_ids
  if (!any(in_chr)) return(rep(0, n_bins * (n_bins - 1L) / 2L))

  d <- dt[in_chr, , drop = FALSE]
  d <- d[d$bin1 != d$bin2, , drop = FALSE]
  if (nrow(d) == 0) return(rep(0, n_bins * (n_bins - 1L) / 2L))

  i <- old2new[d$bin1 + 1L]
  j <- old2new[d$bin2 + 1L]
  sp <- Matrix::sparseMatrix(
    i = pmin(i, j),
    j = pmax(i, j),
    x = d$count,
    dims = c(n_bins, n_bins)
  )
  mat <- as.matrix(sp + Matrix::t(sp))
  mat[upper.tri(mat, diag = FALSE)]
}

write_imputed_offdiag_chr <- function(v, out_file, global_ids) {
  mat <- matrix_long_to_matrix2D_offdiag(as.numeric(v), triangle = "lower")
  mat[!is.finite(mat)] <- 0
  mat <- round(mat)
  keep <- upper.tri(mat, diag = FALSE) & mat > 0
  if (!any(keep)) return(invisible(NULL))

  ij <- which(keep, arr.ind = TRUE)
  out <- data.frame(
    bin1 = global_ids[ij[, 1]],
    bin2 = global_ids[ij[, 2]],
    count = as.integer(mat[keep]),
    stringsAsFactors = FALSE
  )
  write.table(out, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

ref <- read.delim(ref_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
if (!all(c("cell", "depth", "cluster") %in% names(ref))) {
  stop("Reference file must include columns: cell, depth, cluster", call. = FALSE)
}

pair_tbl <- tibble(
  cell_file = list.files(raw_data_dir, full.names = FALSE),
  key = normalize_key(list.files(raw_data_dir, full.names = FALSE))
)

ref2 <- ref %>%
  filter(!is.na(cell), !is.na(cluster), !is.na(depth), depth >= cfg$min_depth) %>%
  mutate(key = normalize_key(cell)) %>%
  left_join(pair_tbl, by = "key") %>%
  filter(!is.na(cell_file)) %>%
  distinct(cell_file, .keep_all = TRUE)

if (nrow(ref2) < 2) stop("Too few matched cells after filtering.", call. = FALSE)

out_root <- path_here("data", "schicimpute_score_compare_from_pairs")
ensure_dir(out_root)

ref_subset_file <- file.path(out_root, "oocyte_zygote_ref_min_depth_5000.tsv")
ref_subset <- ref2 %>%
  transmute(
    cell = cell,
    depth = depth,
    batch = if ("batch" %in% names(ref2)) batch else 1,
    cluster = cluster
  )
write.table(ref_subset, ref_subset_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

imp_dir <- file.path(out_root, "schicimpute_1M")
ensure_dir(imp_dir)

anchors <- read.delim(anchor_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(anchors) <- c("chr", "start", "end", "idx0")
all_chr <- unique(anchors$chr)
if (length(cfg$chr_subset) > 0) all_chr <- intersect(all_chr, cfg$chr_subset)
if (length(all_chr) == 0) stop("No chromosomes selected for imputation.", call. = FALSE)

for (cid in ref2$cell_file) {
  fout <- file.path(imp_dir, cid)
  if (file.exists(fout)) file.remove(fout)
  file.create(fout)
}

for (chr_name in all_chr) {
  cat("\nProcessing ", chr_name, " for scHiCImpute...\n", sep = "")
  chr_ids <- anchors$idx0[anchors$chr == chr_name]
  n_bins <- length(chr_ids)
  if (n_bins < 2) next

  Y_chr <- matrix(0, nrow = n_bins * (n_bins - 1L) / 2L, ncol = nrow(ref2))
  colnames(Y_chr) <- ref2$cell_file

  for (k in seq_len(nrow(ref2))) {
    dt <- read_pair_file(file.path(raw_data_dir, ref2$cell_file[k]))
    if (nrow(dt) > 0) Y_chr[, k] <- build_chr_long_offdiag(dt, chr_ids)
  }

  fit <- HiCImpute::MCMCImpute(
    scHiC = Y_chr,
    bulk = rowSums(Y_chr),
    expected = NULL,
    n = n_bins,
    mc.cores = cfg$mc.cores,
    cutoff = cfg$cutoff,
    niter = cfg$niter,
    burnin = cfg$burnin
  )

  Y_imp_chr <- apply(fit$Impute_SZ, 2, function(x) {
    mat2d <- matrix_long_to_matrix2D_offdiag(x, triangle = "upper")
    matrix2D_to_matrix_long(mat2d, include.diag = FALSE)
  })
  if (is.null(dim(Y_imp_chr))) Y_imp_chr <- matrix(Y_imp_chr, ncol = 1L)
  colnames(Y_imp_chr) <- colnames(Y_chr)

  for (k in seq_len(ncol(Y_imp_chr))) {
    write_imputed_offdiag_chr(Y_imp_chr[, k], file.path(imp_dir, colnames(Y_imp_chr)[k]), chr_ids)
  }

  rm(Y_chr, Y_imp_chr, fit)
  gc(FALSE)
}

schicimpute_scool <- file.path(out_root, "oocyte_zygote_schicimpute_1M.scool")
run_score_cmd(
  score_bin,
  c(
    "cooler",
    "--dset", "oocyte_zygote_schicimpute",
    "--data_dir", imp_dir,
    "--anchor_file", anchor_file,
    "--reference", ref_subset_file,
    "--resolution", "1M",
    "--out", schicimpute_scool
  ),
  "score cooler scHiCImpute"
)

cat("Built scHiCImpute .scool:", schicimpute_scool, "\n")

