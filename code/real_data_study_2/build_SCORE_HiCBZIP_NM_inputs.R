# score_tutorial_1_innerproduct_bhzip_nm_compare.R
# Server-first script: build HiCBZIP-N(M) imputed pairs and .scool for SCORE benchmarking.
# Run: Rscript score_tutorial_1_innerproduct_bhzip_nm_compare.R

required_pkgs <- c("dplyr", "tibble", "jsonlite", "cmdstanr", "parallel")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))

library(dplyr)
library(tibble)
library(jsonlite)
library(cmdstanr)
library(parallel)

if (!exists("path_here", mode = "function")) {
  script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
  script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE)) else getwd()
  source(file.path(script_dir, "..", "_common", "project_paths.R"))
}
source_hicbzip_core(include_score_helpers = TRUE)

score_bin <- Sys.which("score")
if (!nzchar(score_bin)) {
  fallback_candidates <- c(Sys.getenv("HICBZIP_SCORE_BIN", unset = ""))
  existing <- fallback_candidates[file.exists(fallback_candidates)]
  if (length(existing) > 0) score_bin <- existing[[1]]
}
if (!nzchar(score_bin)) stop("`score` command not found.")

raw_data_dir <- path_here("data", "processed", "real_data_study_2", "oocyte_zygote_mm10", "1M")
anchor_file <- path_here("data", "processed", "real_data_study_2", "mm10.genome_split_1M")
ref_file <- path_here("data", "processed", "real_data_study_2", "oocyte_zygote_ref")
if (!dir.exists(raw_data_dir)) stop("Missing raw pair dir: ", raw_data_dir)
if (!file.exists(anchor_file)) stop("Missing anchor file: ", anchor_file)
if (!file.exists(ref_file)) stop("Missing reference file: ", ref_file)

cfg <- list(
  min_depth = 5000,
  nm_r = 1,
  nm_threshold = 2,
  nm_B = 10,
  include_diag = TRUE,
  iter_warmup = 500,
  iter_sampling = 500,
  chains = 1,
  threads_per_chain = 1,
  workers = max(1L, parallel::detectCores() - 2L),
  seed = 123456L,
  force_rebuild_raw_scool = FALSE,
  force_rerun_nm = TRUE,
  chr_subset = character(0)
)

stan_file <- Sys.getenv("HICBZIP_NM_STAN_FILE", unset = path_here("HiCBZIP", "BHZIP_match_normal.stan"))
if (!file.exists(stan_file)) stop("Stan file not found. Set HICBZIP_NM_STAN_FILE. Current: ", stan_file)

# Warning capture (persistent log + snapshot)
log_root <- path_here("results", "real_data_study_2", "input_build_logs")
ensure_dir(log_root)
warning_log_file <- file.path(log_root, "warnings_nm.log")
warning_snapshot_file <- file.path(log_root, "warnings_last50_snapshot.log")
if (file.exists(warning_log_file)) file.remove(warning_log_file)
if (file.exists(warning_snapshot_file)) file.remove(warning_snapshot_file)

log_warning <- function(w) {
  msg <- conditionMessage(w)
  cat(
    format(Sys.time(), "%F %T"),
    " | ", msg, "\n",
    file = warning_log_file, append = TRUE, sep = ""
  )
  invokeRestart("muffleWarning")
}

save_warning_snapshot <- function() {
  w <- warnings()
  if (!is.null(w)) capture.output(w, file = warning_snapshot_file)
}

run_score_cmd <- function(args, label) {
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

normalize_key <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("\\.1m$", "", x, ignore.case = TRUE)
  x
}

read_pair_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(data.frame(bin1 = integer(0), bin2 = integer(0), count = numeric(0)))
  }
  x <- tryCatch(read.table(path, sep = "\t", header = FALSE, stringsAsFactors = FALSE), error = function(e) data.frame())
  if (nrow(x) == 0 || ncol(x) < 3) return(data.frame(bin1 = integer(0), bin2 = integer(0), count = numeric(0)))
  x <- x[, 1:3]
  colnames(x) <- c("bin1", "bin2", "count")
  x$bin1 <- as.integer(x$bin1)
  x$bin2 <- as.integer(x$bin2)
  x$count <- as.numeric(x$count)
  x <- x[!is.na(x$bin1) & !is.na(x$bin2) & !is.na(x$count), , drop = FALSE]
  x
}

write_imputed_chr <- function(v, out_file, global_ids) {
  n_bins <- length(global_ids)
  idx <- 1L
  chunks <- vector("list", n_bins)
  n_chunks <- 0L
  vals <- as.numeric(v)
  vals[!is.finite(vals)] <- 0
  for (i in seq_len(n_bins)) {
    jv <- i:n_bins
    len <- length(jv)
    seg <- vals[idx:(idx + len - 1L)]
    seg[!is.finite(seg)] <- 0
    seg <- round(seg)
    keep <- !is.na(seg) & seg > 0
    if (any(keep)) {
      n_chunks <- n_chunks + 1L
      chunks[[n_chunks]] <- data.frame(
        bin1 = global_ids[i],
        bin2 = global_ids[jv[keep]],
        count = as.integer(seg[keep]),
        stringsAsFactors = FALSE
      )
    }
    idx <- idx + len
  }
  if (n_chunks > 0L) {
    out <- do.call(rbind, chunks[seq_len(n_chunks)])
    write.table(out, file = out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
}

fit_nm_chr <- function(Y_chr, coverage_chr, cfg, mod) {
  N <- nrow(Y_chr)
  K <- ncol(Y_chr)
  n_bins <- as.integer((-1 + sqrt(1 + 8 * N)) / 2)
  if (n_bins * (n_bins + 1) / 2 != N) stop("N pairs mismatch while fitting N(M).")

  list_neighbor_id <- get_neighbor_id(n_bins, r = cfg$nm_r, include_diag = cfg$include_diag)
  lambda_by_row <- rep(mean(coverage_chr), N)

  fit_one_row <- function(i) {
    Y_sim <- Y_chr[i, ]
    if (sum(Y_sim > 0) > cfg$nm_threshold) {
      Y_input <- Y_sim
    } else {
      Y_input <- as.vector(Y_chr[list_neighbor_id[[i]], , drop = FALSE])
    }

    ebe <- get_EBE_ZNB_Gamma_Beta(Y_input, cfg$nm_B, lambda_by_row[i], fix_negative_w = TRUE)
    a_hat <- ebe$a
    b_hat <- ebe$b
    c_hat <- ebe$c
    d_hat <- ebe$d

    a_norm <- log(a_hat / b_hat)
    sigma_mu <- sqrt(1 / a_hat)
    b_norm <- ifelse(c_hat == 0, -1e4, ifelse(d_hat == 0, 1e4, log(c_hat / d_hat)))
    sigma_pi <- ifelse(c_hat == 0 | d_hat == 0, 1e-4,
                       sqrt((c_hat + d_hat)^2 / (c_hat * d_hat * (c_hat + d_hat + 1))))

    data_i <- list(
      N = K,
      Y = Y_sim,
      lambda = coverage_chr,
      a_norm = a_norm,
      sigma_mu = sigma_mu,
      b_norm = b_norm,
      sigma_pi = sigma_pi
    )

    fit <- tryCatch(
      suppressMessages(
        mod$sample(
          data = data_i,
          seed = cfg$seed + i,
          chains = cfg$chains,
          parallel_chains = cfg$chains,
          threads_per_chain = cfg$threads_per_chain,
          iter_warmup = cfg$iter_warmup,
          iter_sampling = cfg$iter_sampling,
          refresh = 0,
          show_messages = FALSE,
          output_dir = tempdir()
        )
      ),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      warning(sprintf("CmdStan fit failed at row %d: %s", i, conditionMessage(fit)))
      return(rep(NA_real_, K))
    }

    mu_tilde <- tryCatch(
      fit$summary(variables = "mu_tilde")$mean,
      error = function(e) numeric(0)
    )

    if (length(mu_tilde) != K) {
      warning(sprintf(
        "mu_tilde length mismatch at row %d: got %d, expected %d; padding with NA",
        i, length(mu_tilde), K
      ))
      out <- rep(NA_real_, K)
      if (length(mu_tilde) > 0) out[seq_len(min(length(mu_tilde), K))] <- as.numeric(mu_tilde[seq_len(min(length(mu_tilde), K))])
      return(out)
    }

    as.numeric(mu_tilde)
  }

  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  res_list <- withCallingHandlers(
    parallel::mclapply(seq_len(N), fit_one_row, mc.cores = cfg$workers, mc.preschedule = TRUE),
    warning = log_warning
  )
  out_mat <- do.call(rbind, lapply(res_list, function(v) as.numeric(v)))
  if (!is.matrix(out_mat)) out_mat <- matrix(as.numeric(out_mat), ncol = K, byrow = TRUE)
  if (ncol(out_mat) != K) {
    warning(sprintf("Final imputation matrix column mismatch: got %d, expected %d; coercing shape", ncol(out_mat), K))
    out_tmp <- matrix(NA_real_, nrow = N, ncol = K)
    rr <- min(nrow(out_mat), N)
    cc <- min(ncol(out_mat), K)
    out_tmp[seq_len(rr), seq_len(cc)] <- out_mat[seq_len(rr), seq_len(cc), drop = FALSE]
    out_mat <- out_tmp
  }
  out_mat
}

set.seed(cfg$seed)
cat("score path:", score_bin, "\n")
cat("stan file:", stan_file, "\n")

ref <- read.delim(ref_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
if (!all(c("cell", "depth", "cluster") %in% names(ref))) stop("Reference must include: cell, depth, cluster")

pair_files <- list.files(raw_data_dir, full.names = FALSE)
pair_files <- pair_files[!grepl("\\.tmp$", pair_files)]
pair_tbl <- tibble(cell_file = pair_files, key = normalize_key(pair_files))

ref2 <- ref %>%
  filter(!is.na(cell), !is.na(cluster), !is.na(depth), depth >= cfg$min_depth) %>%
  mutate(key = normalize_key(cell)) %>%
  left_join(pair_tbl, by = "key") %>%
  filter(!is.na(cell_file)) %>%
  distinct(cell_file, .keep_all = TRUE)
if (nrow(ref2) < 2) stop("Too few matched cells after filtering.")

ref_subset <- ref2 %>%
  transmute(cell = cell,
            depth = depth,
            batch = if ("batch" %in% names(ref2)) batch else 1,
            cluster = cluster)

out_root <- path_here("data", "bhzip_score_compare_from_pairs_nm")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
ref_subset_file <- file.path(out_root, "oocyte_zygote_ref_min_depth_5000.tsv")
write.table(ref_subset, ref_subset_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

raw_scool <- file.path(out_root, "oocyte_zygote_raw_1M.scool")
if (cfg$force_rebuild_raw_scool || !file.exists(raw_scool)) {
  run_score_cmd(c(
    "cooler",
    "--dset", "oocyte_zygote_raw_nm",
    "--data_dir", raw_data_dir,
    "--anchor_file", anchor_file,
    "--reference", ref_subset_file,
    "--resolution", "1M",
    "--out", raw_scool
  ), "score cooler raw (nm workspace)")
}

imp_dir <- file.path(out_root, "bhzip_nm_1M")
dir.create(imp_dir, recursive = TRUE, showWarnings = FALSE)
for (cid in ref2$cell_file) {
  fout <- file.path(imp_dir, cid)
  if (file.exists(fout)) file.remove(fout)
  file.create(fout)
}

anchors <- read.delim(anchor_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(anchors) <- c("chr", "start", "end", "idx0")
all_chr <- unique(anchors$chr)
if (length(cfg$chr_subset) > 0) all_chr <- intersect(all_chr, cfg$chr_subset)
if (length(all_chr) == 0) stop("No chromosomes selected for imputation.")

mod <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))

run_start <- Sys.time()
chr_elapsed <- numeric(0)
total_chr <- length(all_chr)

for (chr_idx in seq_along(all_chr)) {
  chr_name <- all_chr[[chr_idx]]
  chr_start <- Sys.time()
  cat(
    "\n[", chr_idx, "/", total_chr, "] Processing ", chr_name, " for HiCBZIP-N(M)...\n",
    sep = ""
  )
  chr_ids <- anchors$idx0[anchors$chr == chr_name]
  n_bins <- length(chr_ids)
  if (n_bins < 2) {
    cat("Skipped ", chr_name, " (n_bins < 2)\n", sep = "")
    next
  }

  n_pairs <- as.integer(n_bins * (n_bins + 1) / 2)
  cat("n_bins=", n_bins, " | n_pairs=", n_pairs, "\n", sep = "")
  old2new <- integer(max(anchors$idx0) + 1L)
  old2new[chr_ids + 1L] <- seq_len(n_bins)

  Y_chr <- matrix(0, nrow = n_pairs, ncol = nrow(ref2))
  colnames(Y_chr) <- ref2$cell_file

  for (k in seq_len(nrow(ref2))) {
    f <- file.path(raw_data_dir, ref2$cell_file[k])
    dt <- read_pair_file(f)
    if (nrow(dt) == 0) next

    in_chr <- dt$bin1 %in% chr_ids & dt$bin2 %in% chr_ids
    if (!any(in_chr)) next
    d <- dt[in_chr, , drop = FALSE]
    i <- old2new[d$bin1 + 1L]
    j <- old2new[d$bin2 + 1L]
    ii <- pmin(i, j)
    jj <- pmax(i, j)
    rid <- upper_tri_rowmajor_id(ii, jj, n_bins)
    s <- tapply(d$count, rid, sum)
    Y_chr[as.integer(names(s)), k] <- as.numeric(s)
  }

  withCallingHandlers({
    coverage_chr <- compute_bulk_referenced_coverage(Y_chr)$coverage
    Y_imp_chr <- fit_nm_chr(Y_chr, coverage_chr = coverage_chr, cfg = cfg, mod = mod)
    colnames(Y_imp_chr) <- colnames(Y_chr)

    for (k in seq_len(ncol(Y_imp_chr))) {
      write_imputed_chr(
        v = Y_imp_chr[, k],
        out_file = file.path(imp_dir, colnames(Y_imp_chr)[k]),
        global_ids = chr_ids
      )
    }
  }, warning = log_warning)

  rm(Y_chr, Y_imp_chr)
  gc(FALSE)

  chr_sec <- as.numeric(difftime(Sys.time(), chr_start, units = "secs"))
  chr_elapsed <- c(chr_elapsed, chr_sec)
  mean_chr <- mean(chr_elapsed)
  remaining <- total_chr - chr_idx
  eta_sec <- remaining * mean_chr
  cat(
    sprintf(
      "Completed %s in %.1f min | avg/chr %.1f min | ETA %.1f min\n",
      chr_name, chr_sec / 60, mean_chr / 60, eta_sec / 60
    )
  )
}

cat("HiCBZIP-N(M) imputation files written at:", imp_dir, "\n")

bhzip_nm_scool <- file.path(out_root, "oocyte_zygote_bhzip_nm_1M.scool")
run_score_cmd(c(
  "cooler",
  "--dset", "oocyte_zygote_bhzip_nm",
  "--data_dir", imp_dir,
  "--anchor_file", anchor_file,
  "--reference", ref_subset_file,
  "--resolution", "1M",
  "--out", bhzip_nm_scool
), "score cooler bhzip_nm")

cat("Built N(M) scool:", bhzip_nm_scool, "\n")
cat("Raw scool for NM workspace:", raw_scool, "\n")
cat("Reference:", ref_subset_file, "\n")
cat("Total elapsed (min):", round(as.numeric(difftime(Sys.time(), run_start, units = "mins")), 1), "\n")

save_warning_snapshot()
cat("Warning log file:", warning_log_file, "\n")
cat("Warnings snapshot file:", warning_snapshot_file, "\n")
