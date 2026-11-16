# Fit fixed-cohort N(M) from the GB(NB) hyperparameters and whole-cell exposure.
# Run on Linux with CmdStan configured in cmdstanr. Outputs remain unrounded.
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[[1]]))) else getwd()
source(file.path(script_dir, "..", "_common", "project_paths.R"))

out_dir <- Sys.getenv("HICBZIP_SCORE_FIXED150_OUTPUT", path_here("results", "SCORE_oocyte_zygote", "fixed150"))
input_dir <- file.path(out_dir, "gbnb_fits")
metadata_file <- file.path(out_dir, "metadata", "whole_cell_coverage.csv")
stan_file <- path_here("HiCBZIP", "BHZIP_match_normal.stan")
all_chromosomes <- c(paste0("chr", 1:19), "chrX", "chrY")
chromosome_setting <- Sys.getenv("HICBZIP_NM_CHROMOSOMES", "all")
chromosomes <- if (chromosome_setting == "all") all_chromosomes else trimws(strsplit(chromosome_setting, ",", fixed = TRUE)[[1]])
stopifnot(length(chromosomes) > 0L, !anyDuplicated(chromosomes), all(chromosomes %in% all_chromosomes))
require_files(c(metadata_file, stan_file, file.path(input_dir, paste0(chromosomes, ".rds"))),
              label = "fixed-cohort GB(NB) inputs")
cohort <- read.csv(metadata_file, stringsAsFactors = FALSE)
stopifnot(nrow(cohort) == 150L, !anyDuplicated(cohort$cell_file))

cores <- parallel::detectCores()
default_workers <- if (is.finite(cores)) max(1L, cores - 2L) else 1L
workers <- suppressWarnings(as.integer(Sys.getenv("HICBZIP_NM_WORKERS", as.character(default_workers))))
stopifnot(is.finite(workers), workers >= 1L)
cfg <- list(
  chromosomes = chromosomes, prior_mapping = "delta_method",
  iter_warmup = 500L, iter_sampling = 500L, chains = 1L,
  threads_per_chain = 1L, workers = workers, seed = 123456L, block_size = 1000L
)

if (.Platform$OS.type == "windows") stop("Run the fixed-cohort N(M) fit on Linux.")
require_packages("cmdstanr")
output_dir <- ensure_dir(file.path(out_dir, "nm_fits"))
block_dir <- ensure_dir(file.path(output_dir, "blocks"))

# ---- prior_mapping_functions ----
map_gbnb_prior <- function(a, b, c, d) {
  stopifnot(
    length(a) == 1L, length(b) == 1L,
    is.finite(a), is.finite(b), a > 0, b > 0,
    is.finite(c), is.finite(d), c >= 0, d >= 0
  )

  # Delta-method mapping used by the original N(M) implementation.
  mu_mean <- log(a / b)
  mu_sd <- sqrt(1 / a)

  if (c == 0) {
    pi_mean <- -1e4
    pi_sd <- 1e-4
    pi_boundary <- "sampling-zero"
  } else if (d == 0) {
    pi_mean <- 1e4
    pi_sd <- 1e-4
    pi_boundary <- "structural-zero"
  } else {
    pi_mean <- log(c / d)
    pi_sd <- sqrt((c + d)^2 / (c * d * (c + d + 1)))
    pi_boundary <- "interior"
  }

  result <- c(
    mu_mean = mu_mean,
    mu_sd = mu_sd,
    pi_mean = pi_mean,
    pi_sd = pi_sd
  )
  stopifnot(all(is.finite(result)), mu_sd > 0, pi_sd > 0)
  list(values = result, pi_boundary = pi_boundary)
}


# ---- compile_stan ----
Sys.setenv(
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

mod <- cmdstanr::cmdstan_model(
  stan_file,
  cpp_options = list(stan_threads = TRUE)
)


# ---- fitting_functions ----
fit_one_nm_pair <- function(pair_id, input, cfg, mod) {
  hp <- input$hyperparameters[pair_id, , drop = FALSE]
  y <- as.integer(input$Y[pair_id, ])

  # Forked workers inherit the same R tempdir(). Give each pair its own
  # CmdStan output directory so one worker cannot remove another worker's CSV.
  pair_output_dir <- file.path(
    output_dir,
    "_cmdstan_tmp",
    paste0(
      input$chromosome, "_pair_", sprintf("%06d", pair_id),
      "_pid_", Sys.getpid()
    )
  )
  dir.create(pair_output_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(pair_output_dir, recursive = TRUE, force = TRUE), add = TRUE)

  if (hp$c == 0 && hp$d == 0) {
    stop("Invalid Beta boundary at pair ", pair_id)
  }

  prior <- map_gbnb_prior(hp$a, hp$b, hp$c, hp$d)

  if (identical(prior$pi_boundary, "structural-zero")) {
    return(list(
      pair_id = pair_id,
      mean = rep(0, length(y)),
      min_ess_bulk = NA_real_,
      max_mcse_mean = NA_real_,
      prior = prior$values,
      fit_status = "structural-zero-boundary",
      error = NA_character_
    ))
  }

  data_i <- list(
    N = length(y),
    Y = y,
    lambda = as.numeric(input$lambda),
    a_norm = unname(prior$values["mu_mean"]),
    sigma_mu = unname(prior$values["mu_sd"]),
    b_norm = unname(prior$values["pi_mean"]),
    sigma_pi = unname(prior$values["pi_sd"])
  )

  tryCatch({
    fit <- suppressMessages(mod$sample(
      data = data_i,
      seed = cfg$seed + pair_id,
      chains = cfg$chains,
      parallel_chains = cfg$chains,
      threads_per_chain = cfg$threads_per_chain,
      iter_warmup = cfg$iter_warmup,
      iter_sampling = cfg$iter_sampling,
      refresh = 0,
      show_messages = FALSE,
      output_dir = pair_output_dir
    ))

    summary_mu <- fit$summary(variables = "mu_tilde")

    list(
      pair_id = pair_id,
      mean = as.numeric(summary_mu$mean),
      min_ess_bulk = if ("ess_bulk" %in% names(summary_mu)) {
        suppressWarnings(min(summary_mu$ess_bulk, na.rm = TRUE))
      } else {
        NA_real_
      },
      max_mcse_mean = if ("mcse_mean" %in% names(summary_mu)) {
        suppressWarnings(max(summary_mu$mcse_mean, na.rm = TRUE))
      } else {
        NA_real_
      },
      prior = prior$values,
      fit_status = "mcmc",
      error = NA_character_
    )
  }, error = function(e) {
    list(
      pair_id = pair_id,
      mean = rep(NA_real_, length(y)),
      min_ess_bulk = NA_real_,
      max_mcse_mean = NA_real_,
      prior = prior$values,
      fit_status = "error",
      error = conditionMessage(e)
    )
  })
}

fit_nm_chromosome <- function(chr_name, cfg, mod) {
  input_file <- file.path(
    input_dir, paste0(chr_name, ".rds")
  )
  stopifnot(file.exists(input_file))
  input <- readRDS(input_file)
  stopifnot(
    identical(input$chromosome, chr_name),
    nrow(input$Y) == nrow(input$hyperparameters),
    ncol(input$Y) == length(input$lambda),
    all(input$hyperparameters$valid),
    identical(colnames(input$Y), cohort$cell_file),
    isTRUE(all.equal(input$lambda, cohort$lambda_k, tolerance = 1e-14)),
    identical(as.numeric(input$D_k), as.numeric(cohort$D_k)),
    all(is.finite(input$lambda)),
    all(input$lambda > 0),
    abs(sum(input$lambda) - 1) < 1e-12
  )

  final_file <- file.path(output_dir, paste0(chr_name, ".rds"))
  checkpoint_key <- list(
    input_md5 = unname(tools::md5sum(input_file)),
    stan_md5 = unname(tools::md5sum(stan_file)),
    sampler = cfg[c("prior_mapping", "iter_warmup", "iter_sampling",
                    "chains", "threads_per_chain", "seed", "block_size")]
  )
  read_checkpoint <- function(path) {
    result <- readRDS(path)
    if (!identical(result$checkpoint_key, checkpoint_key)) {
      stop("Stale N(M) checkpoint: ", path,
           ". Use a fresh output directory or move aside this chromosome's old fits.")
    }
    result
  }
  if (file.exists(final_file)) return(read_checkpoint(final_file))

  n_pairs <- nrow(input$Y)
  n_cells <- ncol(input$Y)

  if (isTRUE(input$empty_chromosome)) {
    result <- list(
      chromosome = chr_name,
      imputed = matrix(0, n_pairs, n_cells, dimnames = dimnames(input$Y)),
      diagnostics = data.frame(
        pair_id = seq_len(n_pairs), fit_status = "empty-chromosome",
        min_ess_bulk = NA_real_, max_mcse_mean = NA_real_,
        error = NA_character_
      ),
      prior_mapping = cfg$prior_mapping,
      settings = cfg,
      checkpoint_key = checkpoint_key
    )
    saveRDS(result, final_file, compress = "gzip")
    return(result)
  }

  block_starts <- seq.int(1L, n_pairs, by = cfg$block_size)
  for (block_index in seq_along(block_starts)) {
    first <- block_starts[block_index]
    last <- min(n_pairs, first + cfg$block_size - 1L)
    ids <- first:last
    block_file <- file.path(
      block_dir,
      paste0(
        chr_name, "_", cfg$prior_mapping, "_fixed150_pairs_",
        sprintf("%06d", first), "_", sprintf("%06d", last), ".rds"
      )
    )

    if (file.exists(block_file)) {
      stopifnot(identical(read_checkpoint(block_file)$pair_ids, ids))
      cat(chr_name, ": block", first, "to", last, "already exists; skipping.\n")
      next
    }

    cat(chr_name, ": fitting pairs", first, "to", last, "of", n_pairs, "...\n")
    started <- Sys.time()
    fits <- parallel::mclapply(
      ids,
      fit_one_nm_pair,
      input = input,
      cfg = cfg,
      mod = mod,
      mc.cores = cfg$workers,
      mc.preschedule = TRUE
    )

    if (any(!vapply(fits, is.list, logical(1)))) stop("An N(M) worker exited without a fit.")
    errors <- vapply(fits, function(x) !is.na(x$error), logical(1))
    if (any(errors)) {
      failed <- vapply(fits[errors], `[[`, integer(1), "pair_id")
      messages <- vapply(fits[errors], `[[`, character(1), "error")
      stop(
        chr_name, " block failed for pair(s) ",
        paste(failed, collapse = ", "), ": ",
        paste(unique(messages), collapse = " | ")
      )
    }

    block_result <- list(
      chromosome = chr_name,
      pair_ids = ids,
      imputed = do.call(rbind, lapply(fits, `[[`, "mean")),
      diagnostics = data.frame(
        pair_id = vapply(fits, `[[`, integer(1), "pair_id"),
        fit_status = vapply(fits, `[[`, character(1), "fit_status"),
        min_ess_bulk = vapply(fits, `[[`, numeric(1), "min_ess_bulk"),
        max_mcse_mean = vapply(fits, `[[`, numeric(1), "max_mcse_mean"),
        error = vapply(fits, `[[`, character(1), "error"),
        stringsAsFactors = FALSE
      ),
      prior_parameters = do.call(rbind, lapply(fits, `[[`, "prior")),
      elapsed_seconds = as.numeric(
        difftime(Sys.time(), started, units = "secs")
      ),
      prior_mapping = cfg$prior_mapping,
      checkpoint_key = checkpoint_key
    )
    rownames(block_result$imputed) <- as.character(ids)
    colnames(block_result$imputed) <- colnames(input$Y)
    rownames(block_result$prior_parameters) <- as.character(ids)

    stopifnot(
      identical(dim(block_result$imputed), c(length(ids), n_cells)),
      all(is.finite(block_result$imputed)),
      all(block_result$imputed >= 0)
    )
    saveRDS(block_result, block_file, compress = "gzip")
    cat(
      chr_name, ": saved block in",
      round(block_result$elapsed_seconds, 1), "seconds.\n"
    )
  }

  block_files <- file.path(
    block_dir,
    vapply(block_starts, function(first) {
      last <- min(n_pairs, first + cfg$block_size - 1L)
      paste0(
        chr_name, "_", cfg$prior_mapping, "_fixed150_pairs_",
        sprintf("%06d", first), "_", sprintf("%06d", last), ".rds"
      )
    }, character(1))
  )
  stopifnot(all(file.exists(block_files)))
  blocks <- lapply(block_files, read_checkpoint)

  result <- list(
    chromosome = chr_name,
    imputed = do.call(rbind, lapply(blocks, `[[`, "imputed")),
    diagnostics = do.call(rbind, lapply(blocks, `[[`, "diagnostics")),
    prior_parameters = do.call(rbind, lapply(blocks, `[[`, "prior_parameters")),
    prior_mapping = cfg$prior_mapping,
    settings = cfg,
    source_input = input_file,
    checkpoint_key = checkpoint_key
  )

  stopifnot(
    identical(dim(result$imputed), dim(input$Y)),
    all(is.finite(result$imputed)),
    all(result$imputed >= 0),
    identical(result$diagnostics$pair_id, seq_len(n_pairs))
  )
  rownames(result$imputed) <- rownames(input$Y)
  colnames(result$imputed) <- colnames(input$Y)

  saveRDS(result, final_file, compress = "gzip")
  cat(chr_name, ": saved complete fit to", final_file, "\n")
  result
}


# ---- run_selected_chromosomes ----
finite_min_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

finite_max_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

fit_summaries <- vector("list", length(cfg$chromosomes))

for (index in seq_along(cfg$chromosomes)) {
  chr_name <- cfg$chromosomes[index]
  cat("\n=====", chr_name, "=====\n")
  fit <- fit_nm_chromosome(chr_name, cfg, mod)
  fit_summaries[[index]] <- data.frame(
    chromosome = chr_name,
    pairs = nrow(fit$imputed),
    cells = ncol(fit$imputed),
    nonfinite = sum(!is.finite(fit$imputed)),
    nonzero = sum(fit$imputed > 0),
    structural_boundary_pairs = sum(
      fit$diagnostics$fit_status == "structural-zero-boundary"
    ),
    mcmc_pairs = sum(fit$diagnostics$fit_status == "mcmc"),
    minimum_ess_bulk = finite_min_or_na(fit$diagnostics$min_ess_bulk),
    maximum_mcse_mean = finite_max_or_na(fit$diagnostics$max_mcse_mean),
    prior_mapping = fit$prior_mapping,
    stringsAsFactors = FALSE
  )
  rm(fit)
  gc(FALSE)
}

fit_summaries <- do.call(rbind, fit_summaries)
write.csv(
  fit_summaries,
  file.path(
    output_dir,
    "nm_fit_summary.csv"
  ),
  row.names = FALSE
)
fit_summaries
