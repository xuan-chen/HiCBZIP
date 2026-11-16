# Stage manuscript depth-matched inputs from unrounded HiCBZIP fits and
# the preserved HiCImpute output. Math follows the finalized fixed150 scripts.
options(digits = 17)
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
source(file.path(script_dir, "..", "_common", "project_paths.R"))

out_dir <- Sys.getenv("HICBZIP_SCORE_FIXED150_OUTPUT", path_here("results", "SCORE_oocyte_zygote", "fixed150"))
metadata <- file.path(out_dir, "metadata")
coverage <- read.csv(file.path(metadata, "whole_cell_coverage.csv"))
reference <- read.delim(file.path(metadata, "oocyte_zygote_ref_fixed150.tsv"))
cell_key <- function(x) sub("\\.1M$", "", x)
cells <- coverage$cell_file
stopifnot(length(cells) == 150L, !anyDuplicated(cells),
          identical(cell_key(cells), cell_key(reference$cell)),
          all(is.finite(coverage$D_k) & coverage$D_k > 0 & coverage$D_k == floor(coverage$D_k)))
chromosomes <- c(paste0("chr", 1:19), "chrX", "chrY")
models <- commandArgs(TRUE)
if (!length(models)) models <- c("gbnb", "nm", "schicimpute")
stopifnot(!anyDuplicated(models), all(models %in% c("gbnb", "nm", "schicimpute")))
pair_root <- file.path(out_dir, "depth_matched_inputs", "pair_files")

upper_triangle_indices <- function(n_bins) {
  list(
    i = rep(seq_len(n_bins), times = n_bins:1L),
    j = unlist(
      lapply(seq_len(n_bins), function(index) index:n_bins),
      use.names = FALSE
    )
  )
}

# Write one temporary file and rename only after every chromosome succeeds.
# This makes a completed cell file a reliable restart checkpoint.
write_cell_pair_file <- function(values_by_chr, inputs, output_file, integer) {
  temporary_file <- paste0(output_file, ".partial")
  if (file.exists(temporary_file)) unlink(temporary_file)
  if (file.exists(output_file)) unlink(output_file)
  created <- file.create(temporary_file)
  if (!created) stop("Could not create ", temporary_file)

  rows_written <- 0L
  value_sum <- 0
  for (chr_name in names(inputs)) {
    values <- as.numeric(values_by_chr[[chr_name]])
    geometry <- upper_triangle_indices(inputs[[chr_name]]$n_bins)
    stopifnot(length(values) == length(geometry$i))
    keep <- is.finite(values) & values > 0
    if (!any(keep)) next

    output <- data.frame(
      bin1 = inputs[[chr_name]]$anchor_ids[geometry$i[keep]],
      bin2 = inputs[[chr_name]]$anchor_ids[geometry$j[keep]],
      count = if (integer) as.integer(values[keep]) else values[keep]
    )
    write.table(
      output, file = temporary_file, sep = "\t", quote = FALSE,
      row.names = FALSE, col.names = FALSE, append = TRUE
    )
    rows_written <- rows_written + nrow(output)
    value_sum <- value_sum + sum(output$count)
  }

  if (!file.rename(temporary_file, output_file)) {
    stop("Could not finalize ", output_file)
  }
  list(rows = rows_written, total = value_sum)
}

# Allocate a fixed integer total while changing the continuous vector as little
# as possible. Ties are resolved by the fixed chromosome/pair order.
largest_remainder <- function(values, target_total, multiply_first = FALSE) {
  values <- as.numeric(values)
  target_total <- as.integer(round(target_total))
  stopifnot(
    target_total >= 0L, all(is.finite(values)), all(values >= 0),
    sum(values) > 0
  )

  scaled <- if (multiply_first) values * target_total / sum(values) else
    values * (target_total / sum(values))
  allocated <- floor(scaled)
  remainder <- target_total - sum(allocated)
  if (remainder > 0L) {
    fractions <- scaled - allocated
    selected <- order(-fractions, seq_along(fractions))[seq_len(remainder)]
    allocated[selected] <- allocated[selected] + 1
  }
  allocated <- as.integer(allocated)
  stopifnot(sum(allocated) == target_total, all(allocated >= 0L))
  allocated
}

read_saved_pairs <- function(path) {
  require_files(path, "saved HiCImpute output")
  if (file.info(path)$size == 0) {
    return(data.frame(bin1 = integer(), bin2 = integer(), count = numeric()))
  }
  pairs <- read.delim(path, header = FALSE)[, 1:3, drop = FALSE]
  names(pairs) <- c("bin1", "bin2", "count")
  stopifnot(!anyNA(pairs$bin1), !anyNA(pairs$bin2),
            all(is.finite(pairs$count) & pairs$count >= 0))
  pairs[pairs$count > 0, , drop = FALSE]
}

for (model in models) {
  message("Depth matching: ", model)
  continuous_dir <- ensure_dir(file.path(pair_root, "depth_matched_continuous", model))
  integer_dir <- ensure_dir(file.path(pair_root, "depth_matched_integer", model))
  if (model == "schicimpute") {
    saved_dir <- Sys.getenv("HICBZIP_SCORE_HICIMPUTE_PAIRS")
    if (!nzchar(saved_dir)) stop("Set HICBZIP_SCORE_HICIMPUTE_PAIRS to the preserved HiCImpute pair files.")
    require_dirs(saved_dir, "saved HiCImpute pair directory")
    pairs <- lapply(file.path(saved_dir, cells), read_saved_pairs)
    keep <- which(vapply(pairs, nrow, integer(1)) > 0L)
    stopifnot(length(keep) == 146L)
    write.table(reference[keep, , drop = FALSE],
                file.path(metadata, "oocyte_zygote_ref_schicimpute_complete146.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    write.csv(data.frame(cell_file = cells, nonempty_output = seq_along(cells) %in% keep),
              file.path(metadata, "hicimpute_complete_cases.csv"), row.names = FALSE)
  } else {
    fits <- setNames(lapply(file.path(out_dir, "gbnb_fits", paste0(chromosomes, ".rds")), readRDS), chromosomes)
    inputs <- lapply(fits, function(x) x[c("n_bins", "anchor_ids")])
    for (chr in chromosomes) {
      fit <- fits[[chr]]
      stopifnot(identical(chr, fit$chromosome), identical(colnames(fit$Y), cells),
                identical(as.numeric(fit$D_k), as.numeric(coverage$D_k)))
      if (model == "nm") {
        nm <- readRDS(file.path(out_dir, "nm_fits", paste0(chr, ".rds")))
        gb_file <- file.path(out_dir, "gbnb_fits", paste0(chr, ".rds"))
        stopifnot(identical(nm$checkpoint_key$input_md5, unname(tools::md5sum(gb_file))))
        fits[[chr]] <- nm
      }
    }
    matrices <- lapply(fits, "[[", "imputed")
    stopifnot(all(vapply(matrices, function(x) identical(colnames(x), cells), logical(1))))
    lengths <- vapply(matrices, nrow, integer(1))
    chr_index <- rep(chromosomes, times = lengths)
    split_chr <- function(x) split(x, factor(chr_index, levels = chromosomes))
    keep <- seq_along(cells)
  }

  for (k in keep) {
    target <- coverage$D_k[k]
    continuous_file <- file.path(continuous_dir, cells[k])
    integer_file <- file.path(integer_dir, cells[k])
    if (model == "schicimpute") {
      # Keep the saved pair order for integerization ties.
      continuous <- integer <- pairs[[k]]
      continuous$count <- continuous$count * target / sum(continuous$count)
      integer$count <- largest_remainder(integer$count, target, multiply_first = TRUE)
      stopifnot(abs(sum(continuous$count) - target) < 1e-6, sum(integer$count) == target)
      write.table(continuous, continuous_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
      write.table(integer[integer$count > 0, ], integer_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    } else {
      theta <- unlist(lapply(matrices, function(x) x[, k]), use.names = FALSE)
      stopifnot(all(is.finite(theta) & theta >= 0), sum(theta) > 0)
      continuous <- theta * (target / sum(theta))
      integer <- largest_remainder(continuous, target)
      continuous_write <- write_cell_pair_file(split_chr(continuous), inputs, continuous_file, integer = FALSE)
      integer_write <- write_cell_pair_file(split_chr(integer), inputs, integer_file, integer = TRUE)
      stopifnot(abs(continuous_write$total - target) < 1e-6, integer_write$total == target)
    }
  }
  stopifnot(setequal(list.files(continuous_dir), cells[keep]), setequal(list.files(integer_dir), cells[keep]))
}
