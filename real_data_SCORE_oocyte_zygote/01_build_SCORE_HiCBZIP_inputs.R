# Prepare the raw-QC fixed cohort and fit GB(NB) at whole-cell exposure.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))
source_hicbzip_core(include_score_helpers = TRUE)

source_dir <- Sys.getenv("HICBZIP_SCORE_SOURCE", path_here("data", "processed", "SCORE_oocyte_zygote", "source_inputs"))
out_dir <- Sys.getenv("HICBZIP_SCORE_FIXED150_OUTPUT", path_here("results", "SCORE_oocyte_zygote", "fixed150"))
raw_dir <- file.path(source_dir, "oocyte_zygote_mm10", "1M")
anchor_file <- file.path(source_dir, "mm10.genome_split_1M")
reference_file <- file.path(source_dir, "oocyte_zygote_ref")
raw_qc_file <- Sys.getenv(
  "HICBZIP_SCORE_RAW_QC",
  path_here("data", "processed", "SCORE_oocyte_zygote", "source_inputs",
            "oocyte_zygote_raw_pairs_filtered_ref")
)
require_dirs(raw_dir, label = "SCORE raw pair directory")
require_files(c(anchor_file, reference_file, raw_qc_file), label = "SCORE cohort inputs")

all_chromosomes <- c(paste0("chr", 1:19), "chrX", "chrY")
chromosome_setting <- Sys.getenv("HICBZIP_GBNB_CHROMOSOMES", "all")
chromosomes <- if (chromosome_setting == "all") all_chromosomes else trimws(strsplit(chromosome_setting, ",", fixed = TRUE)[[1]])
stopifnot(length(chromosomes) > 0L, !anyDuplicated(chromosomes), all(chromosomes %in% all_chromosomes))

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
cell_key <- function(x) sub("\\.1m$", "", tolower(x))
ref <- read_tsv(reference_file)
qc <- read_tsv(raw_qc_file)
stopifnot(all(c("cell", "depth", "cluster") %in% names(ref)), all(c("cell", "filtered_reason") %in% names(qc)))
raw_files <- list.files(raw_dir, pattern = "\\.1M$", ignore.case = TRUE)
stopifnot(!anyDuplicated(cell_key(raw_files)))
ref <- ref[!is.na(ref$cell) & !is.na(ref$cluster) & !is.na(ref$depth) & ref$depth >= 5000, , drop = FALSE]
ref$cell_file <- raw_files[match(cell_key(ref$cell), cell_key(raw_files))]
ref <- ref[!is.na(ref$cell_file) & !duplicated(ref$cell_file), , drop = FALSE]
keep <- cell_key(qc$cell[which(qc$filtered_reason == "not filtered")])
stopifnot(nrow(ref) == 152L, length(keep) == 150L, !anyDuplicated(keep), all(keep %in% cell_key(ref$cell)))
ref <- ref[cell_key(ref$cell) %in% keep, , drop = FALSE]
cells <- ref$cell_file

anchors <- read.delim(anchor_file, header = FALSE, stringsAsFactors = FALSE)
names(anchors) <- c("chr", "start", "end", "idx0")
stopifnot(!anyNA(anchors$idx0), !anyDuplicated(anchors$idx0), all(all_chromosomes %in% anchors$chr))
anchor_chr <- setNames(anchors$chr, anchors$idx0)
read_intrachromosomal <- function(cell) {
  path <- file.path(raw_dir, cell)
  if (file.info(path)$size == 0) {
    return(data.frame(bin1 = integer(), bin2 = integer(), count = numeric(), chr = character()))
  }
  x <- read.delim(path, header = FALSE, stringsAsFactors = FALSE)
  x <- x[, 1:3, drop = FALSE]
  names(x) <- c("bin1", "bin2", "count")
  left <- unname(anchor_chr[as.character(x$bin1)])
  right <- unname(anchor_chr[as.character(x$bin2)])
  stopifnot(!anyNA(left), !anyNA(right), all(is.finite(x$count) & x$count >= 0 & x$count == floor(x$count)))
  x$chr <- left
  x[left == right & left %in% all_chromosomes, , drop = FALSE]
}
contacts <- setNames(lapply(cells, read_intrachromosomal), cells)
D_k <- vapply(contacts, function(x) sum(x$count), numeric(1))
lambda <- unname(D_k / sum(D_k))
stopifnot(all(D_k > 0), all(is.finite(lambda)), abs(sum(lambda) - 1) < 1e-12)

metadata_dir <- ensure_dir(file.path(out_dir, "metadata"))
write.table(ref[, setdiff(names(ref), "cell_file"), drop = FALSE],
            file.path(metadata_dir, "oocyte_zygote_ref_fixed150.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.csv(data.frame(cell_file = cells, D_k = unname(D_k), lambda_k = lambda),
          file.path(metadata_dir, "whole_cell_coverage.csv"), row.names = FALSE)

if (!identical(Sys.getenv("HICBZIP_SCORE_PREPARE_ONLY"), "1")) {
  input_dir <- ensure_dir(file.path(out_dir, "gbnb_fits"))
  for (chr in chromosomes) {
    message("Fitting GB(NB): ", chr)
    ids <- anchors$idx0[anchors$chr == chr]
    n_bins <- length(ids)
    Y <- matrix(0, n_bins * (n_bins + 1L) / 2L, length(cells),
                dimnames = list(NULL, cells))
    for (k in seq_along(cells)) {
      x <- contacts[[k]]
      x <- x[x$chr == chr, , drop = FALSE]
      if (!nrow(x)) next
      i <- match(x$bin1, ids)
      j <- match(x$bin2, ids)
      sums <- tapply(x$count, upper_tri_rowmajor_id(pmin(i, j), pmax(i, j), n_bins), sum)
      Y[as.integer(names(sums)), k] <- as.numeric(sums)
    }
    empty <- !any(Y > 0)
    if (empty) {
      hp <- data.frame(pair_id = seq_len(nrow(Y)), valid = TRUE, a = NA_real_, b = NA_real_, c = 10, d = 0)
      imputed <- Y
    } else {
      fit <- gbnb_fit_matrix_nested(Y, lambda, n_bins, radius = 1, threshold = 3,
                                    B = 10, include_diagonal = TRUE,
                                    use_global_fallback = TRUE, progress_every = 1000L)
      hp <- fit$hyperparameters
      imputed <- fit$weighted
      stopifnot(all(hp$valid), all(is.finite(imputed)), identical(dim(imputed), dim(Y)))
      rm(fit)
    }
    saveRDS(list(chromosome = chr, n_bins = n_bins, include_diagonal = TRUE,
                 anchor_ids = ids, Y = Y, lambda = lambda, D_k = unname(D_k),
                 hyperparameters = hp, imputed = imputed, empty_chromosome = empty),
            file.path(input_dir, paste0(chr, ".rds")))
  }
}
