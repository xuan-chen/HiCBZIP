## HiCBZIP-GB(NB) empirical Bayes estimation with cell-specific coverage.
##
## Nested-root fit for the mouse oocyte-to-zygote benchmark.
## Returns posterior-weighted imputation and the shared Gamma-Beta parameters.

## ---------------------------------------------------------------------------
## 1. Input checks and zero-count probability
## ---------------------------------------------------------------------------

## Check the count and coverage vectors used by all pair-level functions.
##
## Arguments:
##   Y:      observed counts for one locus pair, one value per cell.
##   lambda: positive cell-specific coverage values in the same cell order.
##
## Returns the numeric vectors in a list after checking their dimensions and
## admissible values. Invalid input stops immediately with an informative error.

gbnb_validate_inputs <- function(Y, lambda) {
  Y <- as.numeric(Y)
  lambda <- as.numeric(lambda)
  if (length(Y) != length(lambda)) {
    stop("Y and lambda must have the same length.")
  }
  if (length(Y) == 0L || any(!is.finite(Y)) || any(Y < 0)) {
    stop("Y must be a nonempty vector of finite, nonnegative counts.")
  }
  if (any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("lambda must contain finite, strictly positive values.")
  }
  list(Y = Y, lambda = lambda)
}

## Negative Binomial probability of a sampling zero after integrating out
## mu ~ Gamma(a, b): P(Y = 0 | nonstructural, lambda) = (b/(b+lambda))^a.
## The log1p form is numerically stable when lambda/b is small.
gbnb_p0 <- function(a, b, lambda) {
  exp(-a * log1p(lambda / b))
}

## Numerically stable calculation of 1 - P(Y = 0 | nonstructural, lambda).
## This quantity is the normalizing probability for the zero-truncated NB.
gbnb_one_minus_p0 <- function(a, b, lambda) {
  -expm1(-a * log1p(lambda / b))
}

## ---------------------------------------------------------------------------
## 2. Cell-specific and pooled zero-truncated moments
## ---------------------------------------------------------------------------

## Calculate the theoretical mean and variance of Y conditional on Y > 0 for
## each supplied lambda_k. The function also returns the untruncated NB zero
## probability p0_k. Here a and b are shared by the locus pair, while lambda_k
## can differ across cells.
gbnb_truncated_moments <- function(a, b, lambda) {
  if (length(a) != 1L || length(b) != 1L ||
      !is.finite(a) || !is.finite(b) || a <= 0 || b <= 0) {
    stop("a and b must be finite positive scalars.")
  }
  lambda <- as.numeric(lambda)
  if (any(!is.finite(lambda)) || any(lambda <= 0)) {
    stop("lambda must contain finite, strictly positive values.")
  }

  x <- lambda / b
  one_minus_p0 <- gbnb_one_minus_p0(a, b, lambda)
  mean_positive <- a * x / one_minus_p0
  second_positive <- (a * x + a * (a + 1) * x^2) / one_minus_p0
  variance_positive <- second_positive - mean_positive^2

  numerical_negative <- variance_positive < 0 &
    variance_positive > -sqrt(.Machine$double.eps) * (1 + abs(second_positive))
  variance_positive[numerical_negative] <- 0

  list(
    mean = mean_positive,
    variance = variance_positive,
    p0 = gbnb_p0(a, b, lambda)
  )
}

## Pool the cell-specific positive-count moments over the cells in which the
## target or borrowed locus pairs have positive observations.
##
## The pooled mean is the average of the cell-specific conditional means. The
## pooled variance is the average within-cell conditional variance plus the
## between-cell variance of the conditional means. At least two positive
## observations are required to define this variance.
gbnb_pooled_moments <- function(a, b, lambda_positive) {
  if (length(lambda_positive) < 2L) {
    stop("At least two positive observations are required for pooled variance.")
  }
  cell_moments <- gbnb_truncated_moments(a, b, lambda_positive)
  m <- cell_moments$mean
  v <- cell_moments$variance
  mbar <- mean(m)
  pooled_variance <- mean(v) + sum((m - mbar)^2) / (length(m) - 1L)
  c(mean = mbar, variance = pooled_variance)
}

## ---------------------------------------------------------------------------
## 3. Nested roots for the pair-level Gamma parameters a_ij and b_ij
## ---------------------------------------------------------------------------

## Inner root: for a fixed Gamma shape a, solve the pooled positive-count mean
## equation for the unique Gamma rate b. The search is performed on log(b) to
## enforce b > 0 and improve numerical stability.
gbnb_solve_b <- function(a, observed_mean, lambda_positive,
                         log_b_range = c(-25, 25), tolerance = 1e-10) {
  if (!is.finite(observed_mean) || observed_mean <= 1) {
    stop("The positive-count mean must exceed one for a finite b root.")
  }
  objective <- function(log_b) {
    unname(gbnb_pooled_moments(
      a = a,
      b = exp(log_b),
      lambda_positive = lambda_positive
    )["mean"]) - observed_mean
  }
  root <- uniroot(
    objective,
    interval = log_b_range,
    extendInt = "downX",
    tol = tolerance
  )$root
  exp(root)
}

## Outer root: estimate a and b from the observed positive counts.
##
## Steps:
##   1. Compute the sample mean and variance among positive counts.
##   2. For each trial a, call gbnb_solve_b() so the mean equation is satisfied.
##   3. Solve the resulting profile-variance equation for a on log(a).
##   4. If the sample variance is outside the attainable profile range, use the
##      appropriate logarithmic-series or Poisson boundary.
##
## Returns a list containing the estimates, observed and fitted moments,
## residuals, profile limits, and a status describing the solution. Fewer than
## three positive counts and an all-ones positive set are returned as invalid
## local fits so the matrix workflow can invoke neighborhood borrowing.
gbnb_fit_ab_nested <- function(Y_positive, lambda_positive,
                               a_range = c(1e-3, 1e3),
                               root_tolerance = 1e-8) {
  Y_positive <- as.numeric(Y_positive)
  lambda_positive <- as.numeric(lambda_positive)
  if (length(Y_positive) != length(lambda_positive)) {
    stop("Y_positive and lambda_positive must have the same length.")
  }
  if (any(!is.finite(Y_positive)) || any(Y_positive <= 0)) {
    stop("Y_positive must contain only finite, strictly positive counts.")
  }
  if (any(!is.finite(lambda_positive)) || any(lambda_positive <= 0)) {
    stop("lambda_positive must contain only finite, strictly positive values.")
  }
  n_positive <- length(Y_positive)
  observed_mean <- if (n_positive > 0L) mean(Y_positive) else NA_real_
  observed_variance <- if (n_positive >= 2L) var(Y_positive) else NA_real_

  if (n_positive < 3L) {
    return(list(
      valid = FALSE,
      a = NA_real_, b = NA_real_, n_positive = n_positive,
      observed_mean = observed_mean, observed_variance = observed_variance,
      status = "too.few.positives"
    ))
  }
  if (observed_mean <= 1) {
    return(list(
      valid = FALSE,
      a = NA_real_, b = NA_real_, n_positive = n_positive,
      observed_mean = observed_mean, observed_variance = observed_variance,
      status = "all.ones"
    ))
  }
  if (length(a_range) != 2L || any(!is.finite(a_range)) ||
      any(a_range <= 0) || a_range[1] >= a_range[2]) {
    stop("a_range must contain two increasing positive values.")
  }

  profile_variance <- function(log_a) {
    # For this candidate a, first enforce the positive-count mean equation.
    a <- exp(log_a)
    b <- gbnb_solve_b(a, observed_mean, lambda_positive)
    unname(gbnb_pooled_moments(a, b, lambda_positive)["variance"])
  }
  profile_residual <- function(log_a) {
    profile_variance(log_a) - observed_variance
  }

  log_a_range <- log(a_range)
  residual_low <- profile_residual(log_a_range[1])
  residual_high <- profile_residual(log_a_range[2])
  residual_tolerance <- root_tolerance * (1 + abs(observed_variance))

  if (!all(is.finite(c(residual_low, residual_high)))) {
    return(list(
      valid = FALSE,
      a = NA_real_, b = NA_real_, n_positive = n_positive,
      observed_mean = observed_mean, observed_variance = observed_variance,
      status = "nonfinite.profile"
    ))
  }
  if (residual_low < residual_high - residual_tolerance) {
    return(list(
      valid = FALSE,
      a = NA_real_, b = NA_real_, n_positive = n_positive,
      observed_mean = observed_mean, observed_variance = observed_variance,
      profile_floor = residual_high + observed_variance,
      profile_ceiling = residual_low + observed_variance,
      status = "nonmonotone.profile"
    ))
  }

  if (abs(residual_low) <= residual_tolerance) {
    # Observed variance is at the high-variance logarithmic-series boundary.
    a <- a_range[1]
    status <- "logseries.limit"
  } else if (abs(residual_high) <= residual_tolerance) {
    # Observed variance is at the low-variance Poisson boundary.
    a <- a_range[2]
    status <- "poisson.limit"
  } else if (residual_low > 0 && residual_high < 0) {
    a <- exp(uniroot(
      profile_residual,
      interval = log_a_range,
      f.lower = residual_low,
      f.upper = residual_high,
      tol = root_tolerance
    )$root)
    status <- "root"
  } else if (residual_high > 0) {
    a <- a_range[2]
    status <- "poisson.limit"
  } else {
    a <- a_range[1]
    status <- "logseries.limit"
  }

  b <- gbnb_solve_b(a, observed_mean, lambda_positive)
  fitted_moments <- gbnb_pooled_moments(a, b, lambda_positive)

  list(
    valid = TRUE,
    a = a,
    b = b,
    n_positive = n_positive,
    observed_mean = observed_mean,
    observed_variance = observed_variance,
    fitted_mean = unname(fitted_moments["mean"]),
    fitted_variance = unname(fitted_moments["variance"]),
    mean_residual = unname(fitted_moments["mean"]) - observed_mean,
    variance_residual = unname(fitted_moments["variance"]) - observed_variance,
    profile_floor = residual_high + observed_variance,
    profile_ceiling = residual_low + observed_variance,
    status = status,
    profile_variance = profile_variance
  )
}

## ---------------------------------------------------------------------------
## 4. Exposure-adjusted fallback without a structural-zero component
## ---------------------------------------------------------------------------

## Refit all counts with a Poisson--Gamma model when the initial zero equation
## would give a negative structural-zero weight. The latent mean is estimated
## as sum(Y_k)/sum(lambda_k), so unequal cell coverage remains in the model.
## Positive extra-Poisson variation gives finite Gamma parameters. Otherwise a
## large-shape Poisson boundary is used. An all-zero input is assigned to the
## structural-zero boundary.
gbnb_fit_poisson_gamma_exposure <- function(Y, lambda, B = 10,
                                            poisson_shape = 1e6,
                                            epsilon = 1e-8) {
  checked <- gbnb_validate_inputs(Y, lambda)
  Y <- checked$Y
  lambda <- checked$lambda
  latent_mean <- sum(Y) / sum(lambda)

  if (!is.finite(latent_mean) || latent_mean <= 0) {
    return(list(
      valid = TRUE,
      a = poisson_shape,
      b = poisson_shape / epsilon,
      w = 1, w_raw = 1, c = B, d = 0,
      pbar0 = 1,
      status = "structural.zero.boundary",
      n_positive = sum(Y > 0),
      latent_mean = 0,
      extra_variance = NA_real_
    ))
  }

  extra_variance <- (
    sum((Y - lambda * latent_mean)^2) - latent_mean * sum(lambda)
  ) / sum(lambda^2)

  if (is.finite(extra_variance) && extra_variance > 0) {
    a <- latent_mean^2 / extra_variance
    b <- a / latent_mean
    status <- "poisson.gamma.refit"
  } else {
    a <- poisson_shape
    b <- poisson_shape / latent_mean
    status <- "poisson.boundary"
  }

  list(
    valid = TRUE,
    a = a, b = b,
    w = 0, w_raw = 0, c = 0, d = B,
    pbar0 = mean(gbnb_p0(a, b, lambda)),
    status = status,
    n_positive = sum(Y > 0),
    latent_mean = latent_mean,
    extra_variance = extra_variance
  )
}

## ---------------------------------------------------------------------------
## 5. Complete fit for one locus pair
## ---------------------------------------------------------------------------

## Fit one pair using the nested-root Gamma estimates and the observed-zero
## equation for the structural-zero weight w.
##
## After estimating a and b from positive counts, the function averages the NB
## sampling-zero probabilities over all cells and solves
##   observed zero fraction = w + (1 - w) * average(p0_k).
## With fixed Beta concentration B, c = B*w and d = B*(1-w). If the raw w is
## materially negative, gbnb_fit_poisson_gamma_exposure() refits the pair with
## w = 0. The returned status records the interior root, boundary, or fallback.
gbnb_fit_pair_nested <- function(Y, lambda, a_range = c(1e-3, 1e3),
                                 B = 10, refit_negative_w = TRUE,
                                 w_tolerance = 1e-10) {
  checked <- gbnb_validate_inputs(Y, lambda)
  Y <- checked$Y
  lambda <- checked$lambda
  positive <- Y > 0

  fit_ab <- gbnb_fit_ab_nested(
    Y_positive = Y[positive],
    lambda_positive = lambda[positive],
    a_range = a_range
  )
  if (!isTRUE(fit_ab$valid)) {
    return(c(fit_ab, list(
      w = NA_real_, w_raw = NA_real_, c = NA_real_, d = NA_real_,
      pbar0 = NA_real_
    )))
  }

  pbar0 <- mean(gbnb_p0(fit_ab$a, fit_ab$b, lambda))
  observed_zero_fraction <- mean(!positive)
  denominator <- 1 - pbar0
  if (!is.finite(denominator) || denominator <= .Machine$double.eps) {
    fit_ab$valid <- FALSE
    fit_ab$w <- NA_real_
    fit_ab$w_raw <- NA_real_
    fit_ab$c <- NA_real_
    fit_ab$d <- NA_real_
    fit_ab$pbar0 <- pbar0
    fit_ab$status <- "zero.equation.failure"
    return(fit_ab)
  }

  w_raw <- (observed_zero_fraction - pbar0) / denominator
  if (isTRUE(refit_negative_w) && is.finite(w_raw) && w_raw < -w_tolerance) {
    fallback <- gbnb_fit_poisson_gamma_exposure(Y, lambda, B = B)
    fallback$nested_status <- fit_ab$status
    fallback$w_raw_before_refit <- w_raw
    return(fallback)
  }

  w <- min(1, max(0, w_raw))
  fit_ab$w <- w
  fit_ab$w_raw <- w_raw
  fit_ab$c <- B * w
  fit_ab$d <- B * (1 - w)
  fit_ab$pbar0 <- pbar0
  fit_ab
}

## ---------------------------------------------------------------------------
## 6. Posterior quantities and the two imputation estimators
## ---------------------------------------------------------------------------

## Given a valid pair-level fit, calculate one row per cell containing:
##   q:              posterior structural-zero probability; q = 0 when Y > 0;
##   posterior_mean: E(mu_ijk | Y_ijk) under the Gamma--Poisson component;
##   theta_weighted: (1 - q_ijk) times the posterior mean; and
##
## When c = 0, the fitted model has no structural-zero component and q = 0.
## When d = 0, it is the structural-zero boundary and q = 1 for observed zeros.
gbnb_impute_pair <- function(Y, lambda, fit) {
  checked <- gbnb_validate_inputs(Y, lambda)
  Y <- checked$Y
  lambda <- checked$lambda
  if (!isTRUE(fit$valid) || any(!is.finite(c(fit$a, fit$b, fit$c, fit$d)))) {
    stop("fit must contain finite, valid a, b, c, and d values.")
  }

  if (fit$c <= 0) {
    q_zero <- rep(0, length(Y))
  } else if (fit$d <= 0) {
    q_zero <- rep(1, length(Y))
  } else {
    p0 <- gbnb_p0(fit$a, fit$b, lambda)
    q_zero <- fit$c / (fit$c + fit$d * p0)
  }
  q <- ifelse(Y == 0, q_zero, 0)
  posterior_mean <- (fit$a + Y) / (fit$b + lambda)

  data.frame(
    Y = Y,
    lambda = lambda,
    q = q,
    posterior_mean = posterior_mean,
    theta_weighted = (1 - q) * posterior_mean
  )
}

## ---------------------------------------------------------------------------
## 7. Triangular locus-pair indexing and neighborhood construction
## ---------------------------------------------------------------------------

## Recover the number of genomic bins from the number of stored triangular
## locus pairs. The formula differs depending on whether diagonal pairs are
## included. An incompatible row count triggers an error.
gbnb_infer_n_bins <- function(n_pairs, include_diagonal = TRUE) {
  if (include_diagonal) {
    n <- (-1 + sqrt(1 + 8 * n_pairs)) / 2
    expected <- function(x) x * (x + 1) / 2
  } else {
    n <- (1 + sqrt(1 + 8 * n_pairs)) / 2
    expected <- function(x) x * (x - 1) / 2
  }
  n <- as.integer(round(n))
  if (expected(n) != n_pairs) {
    stop("n_pairs is incompatible with the requested triangular format.")
  }
  n
}

## Build a list of valid neighboring locus-pair row IDs for every (i,j) pair.
## Neighbors satisfy Manhattan distance |delta_i| + |delta_j| <= radius and
## remain inside the stored upper triangle. The target pair is included.
gbnb_neighbor_ids <- function(n_bins, radius = 1, include_diagonal = TRUE) {
  coordinates <- do.call(rbind, lapply(seq_len(n_bins), function(i) {
    first_j <- if (include_diagonal) i else i + 1L
    if (first_j > n_bins) return(NULL)
    cbind(i = i, j = seq.int(first_j, n_bins))
  }))
  coordinates <- as.data.frame(coordinates)
  coordinates$id <- seq_len(nrow(coordinates))
  id_map <- matrix(NA_integer_, nrow = n_bins, ncol = n_bins)
  id_map[cbind(coordinates$i, coordinates$j)] <- coordinates$id

  offsets <- expand.grid(di = -radius:radius, dj = -radius:radius)
  offsets <- offsets[abs(offsets$di) + abs(offsets$dj) <= radius, , drop = FALSE]

  lapply(seq_len(nrow(coordinates)), function(index) {
    neighbor_i <- coordinates$i[index] + offsets$di
    neighbor_j <- coordinates$j[index] + offsets$dj
    valid <- neighbor_i >= 1 & neighbor_i <= n_bins &
      neighbor_j >= 1 & neighbor_j <= n_bins
    ids <- id_map[cbind(neighbor_i[valid], neighbor_j[valid])]
    ids[!is.na(ids)]
  })
}

## ---------------------------------------------------------------------------
## 8. Assemble the data used to fit one target pair
## ---------------------------------------------------------------------------

## Use the target row alone when it has at least 'threshold' positive cells.
## Otherwise pool its precomputed neighboring rows. Because R flattens matrices
## column by column, rep(lambda, each = number of selected rows) aligns every
## pooled count with the coverage of its original cell.
##
## Returns the flattened count and coverage vectors, selected pair IDs, number
## of pooled rows, and whether the fit uses the target or its neighborhood.
gbnb_prepare_fit_data <- function(Y_matrix, lambda, pair_id, neighbor_ids,
                                  threshold = 3) {
  Y_matrix <- as.matrix(Y_matrix)
  lambda <- as.numeric(lambda)
  if (length(lambda) != ncol(Y_matrix)) {
    stop("lambda must have one value for each cell column.")
  }
  target <- as.numeric(Y_matrix[pair_id, ])
  use_neighborhood <- sum(target > 0) < threshold
  fit_pair_ids <- if (use_neighborhood) neighbor_ids[[pair_id]] else pair_id
  block <- Y_matrix[fit_pair_ids, , drop = FALSE]

  list(
    Y = as.vector(block),
    lambda = rep(lambda, each = nrow(block)),
    pair_ids = fit_pair_ids,
    n_neighbor_pairs = nrow(block),
    source = if (use_neighborhood) "neighborhood" else "target"
  )
}

## Create the limiting fit for a target pair that is zero in every cell after
## local fitting has failed. This records w = 1, c = B, and d = 0, so observed
## zeros receive q = 1 and weighted imputed values of zero.
gbnb_structural_boundary_fit <- function(B = 10, poisson_shape = 1e6,
                                         epsilon = 1e-8) {
  list(
    valid = TRUE,
    a = poisson_shape,
    b = poisson_shape / epsilon,
    w = 1, w_raw = 1, c = B, d = 0,
    pbar0 = 1,
    n_positive = 0,
    status = "structural.zero.boundary"
  )
}

## ---------------------------------------------------------------------------
## 9. Fit selected rows of a complete chromosome matrix
## ---------------------------------------------------------------------------

## Apply the pair-level workflow to specified locus-pair rows of a chromosome.
## If pair_ids is NULL, all rows are fitted. For each selected target, the
## function:
##   1. chooses target-only or neighborhood data;
##   2. estimates pair-level hyperparameters with cell-specific lambda_k;
##   3. assigns an unresolved all-zero target to the structural boundary;
##   4. optionally uses one chromosome-wide fit for another unresolved pair; and
##   5. stores weighted, q, posterior-mean, and diagnostic outputs.
##
## n_bins and include_diagonal must describe the geometry used to construct
## Y_matrix. Their implied triangular row count is checked before fitting. This
## explicit check is necessary because the same triangular number can describe
## n bins with the diagonal or n + 1 bins without the diagonal.
##
## Returned matrices have dimensions length(pair_ids) by number of cells. The
## hyperparameters data frame has one row per selected target pair and records
## the data source, mathematical status, sparse-data counts, estimates, and
## moment residuals. The finalized mouse benchmark enables the chromosome-wide
## fallback with use_global_fallback = TRUE.
gbnb_fit_matrix_nested <- function(Y_matrix, lambda, n_bins, pair_ids = NULL,
                                   radius = 1, threshold = 3, B = 10,
                                   include_diagonal = TRUE,
                                   a_range = c(1e-3, 1e3),
                                   use_global_fallback = FALSE,
                                   verbose = TRUE,
                                   progress_every = 100L) {
  Y_matrix <- as.matrix(Y_matrix)
  storage.mode(Y_matrix) <- "double"
  lambda <- as.numeric(lambda)
  if (length(lambda) != ncol(Y_matrix) || any(!is.finite(lambda)) ||
      any(lambda <= 0)) {
    stop("lambda must contain one finite positive value per cell.")
  }
  if (length(n_bins) != 1L || !is.finite(n_bins) || n_bins < 2 ||
      n_bins != as.integer(n_bins)) {
    stop("n_bins must be one integer greater than or equal to 2.")
  }
  if (length(progress_every) != 1L || !is.finite(progress_every) ||
      progress_every < 1 || progress_every != as.integer(progress_every)) {
    stop("progress_every must be one positive integer.")
  }
  progress_every <- as.integer(progress_every)
  n_bins <- as.integer(n_bins)
  expected_n_pairs <- if (include_diagonal) {
    n_bins * (n_bins + 1L) / 2L
  } else {
    n_bins * (n_bins - 1L) / 2L
  }
  if (nrow(Y_matrix) != expected_n_pairs) {
    stop(
      "Y_matrix has ", nrow(Y_matrix), " rows, but n_bins = ", n_bins,
      " and include_diagonal = ", include_diagonal, " require ",
      expected_n_pairs, " triangular locus-pair rows."
    )
  }
  if (is.null(pair_ids)) pair_ids <- seq_len(nrow(Y_matrix))
  pair_ids <- as.integer(pair_ids)
  if (any(pair_ids < 1L | pair_ids > nrow(Y_matrix))) {
    stop("pair_ids contains an out-of-range value.")
  }

  neighbors <- gbnb_neighbor_ids(n_bins, radius, include_diagonal)
  n_selected <- length(pair_ids)
  n_cells <- ncol(Y_matrix)

  weighted <- matrix(NA_real_, n_selected, n_cells)
  q <- matrix(NA_real_, n_selected, n_cells)
  posterior_mean <- matrix(NA_real_, n_selected, n_cells)
  hyper <- vector("list", n_selected)
  global_fit <- NULL

  for (position in seq_along(pair_ids)) {
    pair_id <- pair_ids[position]
    fit_data <- gbnb_prepare_fit_data(
      Y_matrix, lambda, pair_id, neighbors, threshold
    )
    fit <- gbnb_fit_pair_nested(
      fit_data$Y,
      fit_data$lambda,
      a_range = a_range,
      B = B
    )
    fit_source <- fit_data$source

    ## Preserve the explicit all-zero boundary before considering a global fit.
    ## Otherwise the pooled chromosome fit would replace the intended q = 1
    ## result for a target pair containing no observed contacts.
    if (!isTRUE(fit$valid) && sum(Y_matrix[pair_id, ] > 0) == 0L) {
      fit <- gbnb_structural_boundary_fit(B = B)
      fit_source <- "structural-boundary"
    }

    ## The legacy GB(NB) workflow used one chromosome-wide estimate when the
    ## target or its radius-1 neighborhood did not yield a usable fit. Reuse
    ## that rule here with the corrected unequal-coverage estimator. Only the
    ## hyperparameters are borrowed: imputation below still uses the target
    ## pair's own counts and the original cell-specific lambda values.
    if (!isTRUE(fit$valid) && use_global_fallback) {
      if (is.null(global_fit)) {
        global_fit <- gbnb_fit_pair_nested(
          as.vector(Y_matrix),
          rep(lambda, each = nrow(Y_matrix)),
          a_range = a_range,
          B = B
        )
      }
      if (isTRUE(global_fit$valid)) {
        fit <- global_fit
        fit_source <- "global"
      }
    }

    if (isTRUE(fit$valid)) {
      imputation <- gbnb_impute_pair(Y_matrix[pair_id, ], lambda, fit)
      weighted[position, ] <- imputation$theta_weighted
      q[position, ] <- imputation$q
      posterior_mean[position, ] <- imputation$posterior_mean
    }

    hyper[[position]] <- data.frame(
      pair_id = pair_id,
      valid = isTRUE(fit$valid),
      fit_source = fit_source,
      status = if (!is.null(fit$status)) fit$status else NA_character_,
      n_positive_target = sum(Y_matrix[pair_id, ] > 0),
      n_positive_fit = if (!is.null(fit$n_positive)) fit$n_positive else NA_integer_,
      n_neighbor_pairs = fit_data$n_neighbor_pairs,
      a = if (!is.null(fit$a)) fit$a else NA_real_,
      b = if (!is.null(fit$b)) fit$b else NA_real_,
      w = if (!is.null(fit$w)) fit$w else NA_real_,
      w_raw = if (!is.null(fit$w_raw)) fit$w_raw else NA_real_,
      c = if (!is.null(fit$c)) fit$c else NA_real_,
      d = if (!is.null(fit$d)) fit$d else NA_real_,
      mean_residual = if (!is.null(fit$mean_residual)) fit$mean_residual else NA_real_,
      variance_residual = if (!is.null(fit$variance_residual)) {
        fit$variance_residual
      } else {
        NA_real_
      }
    )

    if (verbose &&
        (position %% progress_every == 0L || position == n_selected)) {
      message("Completed ", position, " of ", n_selected, " selected pairs.")
    }
  }

  rownames(weighted) <- rownames(q) <-
    rownames(posterior_mean) <- as.character(pair_ids)
  colnames(weighted) <- colnames(q) <-
    colnames(posterior_mean) <- colnames(Y_matrix)

  list(
    pair_ids = pair_ids,
    weighted = weighted,
    q = q,
    posterior_mean = posterior_mean,
    hyperparameters = do.call(rbind, hyper),
    settings = list(
      n_bins = n_bins,
      radius = radius,
      threshold = threshold,
      B = B,
      include_diagonal = include_diagonal,
      a_range = a_range,
      use_global_fallback = use_global_fallback,
      progress_every = progress_every
    )
  )
}
