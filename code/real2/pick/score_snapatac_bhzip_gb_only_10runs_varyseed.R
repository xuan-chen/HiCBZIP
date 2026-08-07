# score_snapatac_bhzip_gb_only_10runs_varyseed.R
# Run with: Rscript score_snapatac_bhzip_gb_only_10runs_varyseed.R
# Purpose: launch SnapATAC on BHZIP-GB only (no raw), 10 seeds.

score_bin <- Sys.which("score")
if (!nzchar(score_bin)) {
  fallback <- c("/home/xchen/miniconda3/envs/score/bin/score", "C:/Users/67402/.conda/envs/score_py39/Scripts/score.exe")
  hit <- fallback[file.exists(fallback)]
  if (length(hit) > 0) score_bin <- hit[[1]]
}
if (!nzchar(score_bin)) stop("`score` command not found.")

root_dir <- file.path("data", "bhzip_score_compare_from_pairs")
bhzip_scool <- file.path(root_dir, "oocyte_zygote_bhzip_1M.scool")
ref_file <- file.path(root_dir, "oocyte_zygote_ref_min_depth_5000.tsv")
if (!file.exists(bhzip_scool)) stop("Missing BHZIP-GB scool: ", bhzip_scool)
if (!file.exists(ref_file)) stop("Missing reference: ", ref_file)

n_runs <- 10L
seed_start <- 2026L
embedding_alg <- "SnapATAC"

# Toggle settings:
# FALSE => keep IDF (standard SnapATAC behavior)
# TRUE  => add --snapatac_no_idf
snapatac_no_idf <- TRUE
use_xy <- TRUE
no_viz <- FALSE

# Change this tag to isolate output namespace.
run_tag <- if (snapatac_no_idf) "snapatac_no_idf_gb_only" else "snapatac_idf_gb_only"

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
  }

  list(status = as.integer(status), stdout = out, stderr = err)
}

run_one <- function(run_idx, seed) {
  dset_name <- sprintf("oocyte_zygote_bhzip_pairs_%s_r%02d", run_tag, run_idx)
  args <- c(
    "embed",
    "--dset", dset_name,
    "--resolution", "1M",
    "--scool", bhzip_scool,
    "--reference", ref_file,
    "--embedding_algs", embedding_alg,
    "--min_depth", "5000",
    "--seed", as.character(seed)
  )

  if (isTRUE(use_xy)) args <- c(args, "--use_xy")
  if (isTRUE(no_viz)) args <- c(args, "--no_viz")
  if (isTRUE(snapatac_no_idf)) args <- c(args, "--snapatac_no_idf")

  cmd <- run_score_cmd(score_bin, args, label = paste("score embed bhzip_gb run", run_idx))

  fail_msg <- NA_character_
  if (!identical(cmd$status, 0L)) {
    msg <- c(cmd$stderr, cmd$stdout)
    msg <- msg[nzchar(msg)]
    fail_msg <- if (length(msg) > 0) tail(msg, 1) else sprintf("Non-zero return code: %d", cmd$status)
  }

  list(
    run = run_idx,
    seed = seed,
    dset = dset_name,
    status = cmd$status,
    success = identical(cmd$status, 0L),
    failure_msg = fail_msg
  )
}

all_status <- vector("list", n_runs)
for (i in seq_len(n_runs)) {
  seed_i <- seed_start + i - 1L
  all_status[[i]] <- run_one(i, seed_i)
}

status_df <- dplyr::bind_rows(all_status)

cat("\n================ RUN SUMMARY ================\n")
print(status_df)

failed <- dplyr::filter(status_df, !success)
if (nrow(failed) > 0) {
  cat("\nFailed runs:\n")
  print(failed)
  warning("Some bhzip_gb SnapATAC runs failed. See log output above.")
} else {
  cat("\nAll bhzip_gb SnapATAC runs completed successfully.\n")
}
