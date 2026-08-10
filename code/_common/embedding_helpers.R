get_chr_groups <- function(feature_names) {
  if (is.null(feature_names) || !any(grepl("::", feature_names, fixed = TRUE))) {
    return(list(all = seq_along(feature_names)))
  }
  chr <- sub("::.*$", "", feature_names)
  split(seq_along(feature_names), chr)
}

infer_n_from_pairs <- function(n_pairs) {
  n <- (-1 + sqrt(1 + 8 * n_pairs)) / 2
  if (!isTRUE(all.equal(n, round(n)))) stop("Invalid pair length for triangular matrix: ", n_pairs)
  as.integer(round(n))
}

vec_upper_rowmajor_to_mat <- function(v) {
  n <- infer_n_from_pairs(length(v))
  M <- matrix(0, n, n)
  idx <- 1L
  for (i in seq_len(n)) {
    len <- n - i + 1L
    jseq <- i:n
    M[i, jseq] <- v[idx:(idx + len - 1L)]
    idx <- idx + len
  }
  M[lower.tri(M)] <- t(M)[lower.tri(M)]
  M
}

mat_to_upper_rowmajor <- function(M) {
  n <- nrow(M)
  out <- numeric(n * (n + 1) / 2)
  idx <- 1L
  for (i in seq_len(n)) {
    len <- n - i + 1L
    out[idx:(idx + len - 1L)] <- M[i, i:n]
    idx <- idx + len
  }
  out
}

box_smooth_2d <- function(M, radius = 1) {
  if (radius <= 0) return(M)
  n <- nrow(M)
  out <- matrix(0, n, n)
  for (i in seq_len(n)) {
    i1 <- max(1, i - radius)
    i2 <- min(n, i + radius)
    for (j in seq_len(n)) {
      j1 <- max(1, j - radius)
      j2 <- min(n, j + radius)
      out[i, j] <- mean(M[i1:i2, j1:j2])
    }
  }
  out
}

prep_cell_feature_matrix <- function(Y_feature_by_cell) {
  X <- t(as.matrix(Y_feature_by_cell))
  rs <- rowSums(X)
  rs[rs == 0] <- 1
  X <- sweep(X, 1, rs, "/")
  X <- log1p(X * 1e6)
  sds <- apply(X, 2, sd)
  X <- X[, sds > 0, drop = FALSE]
  X
}

run_two_stage_pca <- function(Y_feature_by_cell, feature_names, npc = 20, pcs_per_chr = 5, max_features_per_chr = 8000) {
  if (length(max_features_per_chr) == 0 || is.null(max_features_per_chr) ||
      !is.finite(max_features_per_chr) || max_features_per_chr < 2) {
    max_features_per_chr <- Inf
  } else {
    max_features_per_chr <- as.integer(max_features_per_chr[1])
  }

  chr_groups <- get_chr_groups(feature_names)

  chr_pcs <- lapply(names(chr_groups), function(chr) {
    idx <- chr_groups[[chr]]
    idx <- idx[idx <= nrow(Y_feature_by_cell)]
    if (length(idx) < 2) return(NULL)

    Yc <- Y_feature_by_cell[idx, , drop = FALSE]

    if (is.finite(max_features_per_chr) && nrow(Yc) > max_features_per_chr) {
      s <- rowSums(Yc)
      keep <- order(s, decreasing = TRUE)[seq_len(max_features_per_chr)]
      Yc <- Yc[keep, , drop = FALSE]
    }

    Xc <- t(as.matrix(Yc))
    rm(Yc); gc(FALSE)

    rs <- rowSums(Xc)
    rs[rs == 0] <- 1
    Xc <- sweep(Xc, 1, rs, "/")
    Xc <- log1p(Xc * 1e6)

    sdc <- apply(Xc, 2, sd)
    Xc <- Xc[, sdc > 0, drop = FALSE]
    if (ncol(Xc) < 2) return(NULL)

    pca_chr <- prcomp(Xc, center = TRUE, scale. = TRUE)
    npc_chr <- min(pcs_per_chr, ncol(pca_chr$x))
    out <- pca_chr$x[, seq_len(npc_chr), drop = FALSE]
    rm(Xc, pca_chr); gc(FALSE)
    out
  })

  chr_pcs <- chr_pcs[!vapply(chr_pcs, is.null, logical(1))]
  if (length(chr_pcs) == 0) stop("Two-stage PCA failed: no chromosome-level PCs.")

  X_genome <- do.call(cbind, chr_pcs)
  pca_genome <- prcomp(X_genome, center = TRUE, scale. = TRUE)
  npc_i <- min(npc, ncol(pca_genome$x))
  list(
    pcs_for_cluster = pca_genome$x[, seq_len(npc_i), drop = FALSE],
    emb2d = pca_genome$x[, 1:2, drop = FALSE]
  )
}

schicluster_like_embed <- function(Y_feature_by_cell, cfg, force_pipeline = "pca_umap") {
  feat_names <- rownames(Y_feature_by_cell)
  chr_groups <- get_chr_groups(feat_names)
  if (length(chr_groups) == 1 && names(chr_groups)[1] == "all") {
    warning("No chromosome-prefixed features (chr::pair). Falling back to pca_umap.")
    return(NULL)
  }

  K <- ncol(Y_feature_by_cell)
  chr_pcs <- list()

  for (chr in names(chr_groups)) {
    idx_chr <- chr_groups[[chr]]
    Yc <- Y_feature_by_cell[idx_chr, , drop = FALSE]
    n_pairs <- nrow(Yc)

    ok_chr <- TRUE
    n_bins <- tryCatch(infer_n_from_pairs(n_pairs), error = function(e) { ok_chr <<- FALSE; NA_integer_ })
    if (!ok_chr || is.na(n_bins) || n_bins < 2) next

    X_chr <- matrix(0, nrow = K, ncol = n_pairs)

    for (k in seq_len(K)) {
      v <- as.numeric(Yc[, k])
      M <- vec_upper_rowmajor_to_mat(v)
      rs <- rowSums(M)
      rs[rs <= 0] <- 1
      d <- sqrt(rs)
      M <- M / outer(d, d)
      M[!is.finite(M)] <- 0

      M_box <- box_smooth_2d(M, radius = cfg$schi_box_radius)
      pr <- rowSums(M_box)
      pr[pr <= 0] <- 1
      P <- M_box / pr
      M_rw <- (1 - cfg$schi_rw_alpha) * M_box + cfg$schi_rw_alpha * (P %*% M_box)

      thr <- as.numeric(stats::quantile(M_rw, probs = cfg$schi_top_quantile, na.rm = TRUE))
      M_rw[M_rw < thr] <- 0

      X_chr[k, ] <- mat_to_upper_rowmajor(M_rw)
    }

    sds <- apply(X_chr, 2, sd)
    X_chr <- X_chr[, sds > 0, drop = FALSE]
    if (ncol(X_chr) < 2) next

    pca_chr <- prcomp(X_chr, center = TRUE, scale. = TRUE)
    npc_chr <- min(cfg$schi_pcs_per_chr, ncol(pca_chr$x))
    chr_pcs[[chr]] <- pca_chr$x[, seq_len(npc_chr), drop = FALSE]
  }

  if (length(chr_pcs) == 0) {
    warning("No chromosome produced valid PCs in schicluster_like. Falling back to pca_umap.")
    return(NULL)
  }

  X_genome <- do.call(cbind, chr_pcs)
  pca_genome <- prcomp(X_genome, center = TRUE, scale. = TRUE)
  npc_i <- min(cfg$npc, ncol(pca_genome$x))
  pcs <- pca_genome$x[, seq_len(npc_i), drop = FALSE]

  if (identical(force_pipeline, "pca_umap")) {
    uwot::umap(
      pcs,
      n_neighbors = cfg$umap_n_neighbors,
      min_dist = cfg$umap_min_dist,
      metric = cfg$umap_metric,
      n_components = 2,
      verbose = FALSE,
      ret_model = FALSE
    )
  } else {
    pcs[, 1:2, drop = FALSE]
  }
}

robust_clip_pcs <- function(pcs, clip_z = 5) {
  if (is.null(clip_z) || !is.finite(clip_z) || clip_z <= 0) return(pcs)
  s <- apply(pcs, 2, stats::mad, constant = 1, na.rm = TRUE)
  s[!is.finite(s) | s <= 0] <- 1
  z <- sweep(pcs, 2, s, "/")
  pmax(pmin(z, clip_z), -clip_z)
}


build_cis_visibility_matrix <- function(Y_feature_by_cell, feature_names, max_bins = Inf) {
  chr_groups <- get_chr_groups(feature_names)
  if (length(chr_groups) == 1 && names(chr_groups)[1] == "all") {
    stop("1d_pca requires chromosome-tagged features like chr::pair")
  }

  K <- ncol(Y_feature_by_cell)
  vis_blocks <- lapply(names(chr_groups), function(chr) {
    idx <- chr_groups[[chr]]
    idx <- idx[idx <= nrow(Y_feature_by_cell)]
    if (length(idx) < 2) return(NULL)

    Yc <- Y_feature_by_cell[idx, , drop = FALSE]
    ok <- TRUE
    n_bins <- tryCatch(infer_n_from_pairs(nrow(Yc)), error = function(e) { ok <<- FALSE; NA_integer_ })
    if (!ok || is.na(n_bins) || n_bins < 2) return(NULL)

    out <- matrix(0, nrow = K, ncol = n_bins)
    for (k in seq_len(K)) {
      v <- as.numeric(Yc[, k])
      M <- vec_upper_rowmajor_to_mat(v)
      out[k, ] <- rowSums(M)
    }
    colnames(out) <- paste0(chr, ":bin", seq_len(n_bins))
    out
  })

  vis_blocks <- vis_blocks[!vapply(vis_blocks, is.null, logical(1))]
  if (length(vis_blocks) == 0) stop("1d_pca failed: no valid chromosome blocks")

  X <- do.call(cbind, vis_blocks)
  if (is.finite(max_bins) && ncol(X) > max_bins) {
    s <- colSums(X)
    keep <- order(s, decreasing = TRUE)[seq_len(as.integer(max_bins))]
    X <- X[, keep, drop = FALSE]
  }

  X <- log1p(X)
  sds <- apply(X, 2, sd)
  X <- X[, sds > 0, drop = FALSE]
  X
}

build_innerproduct_distance_paper <- function(Y_feature_by_cell, feature_names, n_strata = 10L) {
  chr_groups <- get_chr_groups(feature_names)
  if (length(chr_groups) == 1 && names(chr_groups)[1] == "all") {
    stop("paper_innerproduct_mds requires chromosome-tagged features like chr::pair")
  }

  if (is.null(n_strata) || length(n_strata) == 0 || !is.finite(n_strata[1]) || n_strata[1] < 1) {
    n_strata <- 10L
  } else {
    n_strata <- as.integer(n_strata[1])
  }

  n_cells <- ncol(Y_feature_by_cell)
  if (is.null(n_cells) || !is.finite(n_cells) || n_cells < 2) {
    stop("paper_innerproduct_mds requires at least 2 cells; got ", n_cells)
  }

  d_blocks <- lapply(names(chr_groups), function(chr) {
    idx <- chr_groups[[chr]]
    idx <- idx[idx <= nrow(Y_feature_by_cell)]
    if (length(idx) < 3) return(NULL)

    Yc <- Y_feature_by_cell[idx, , drop = FALSE]
    ok <- TRUE
    n_bins <- tryCatch(infer_n_from_pairs(nrow(Yc)), error = function(e) { ok <<- FALSE; NA_integer_ })
    if (!ok || is.na(n_bins) || n_bins < 2) return(NULL)

    Xc <- as.matrix(Matrix::t(Yc))
    if (!is.numeric(Xc)) storage.mode(Xc) <- "double"
    n_use <- min(n_strata, n_bins)
    if (n_use < 1) return(NULL)

    strata_feats <- lapply(seq_len(n_use), function(d) {
      j <- d - 1L
      m <- n_bins - j
      cols <- integer(m)
      for (k in seq_len(m)) {
        i1 <- k
        j1 <- k + j
        cols[k] <- as.integer((i1 - 1L) * (2L * n_bins - i1 + 2L) / 2L + (j1 - i1 + 1L))
      }
      S <- Xc[, cols, drop = FALSE]
      mu <- rowMeans(S)
      sdv <- apply(S, 1, stats::sd)
      sdv[!is.finite(sdv) | sdv <= 0] <- 1
      Z <- (S - mu) / sdv
      Z[!is.finite(Z)] <- 0
      Z
    })

    Zcat <- do.call(cbind, strata_feats)
    if (is.null(Zcat) || ncol(Zcat) < 1) return(NULL)

    sim <- (Zcat %*% t(Zcat)) / ncol(Zcat)
    sim <- as.matrix(sim)
    sim[sim > 1] <- 1
    sim[sim < -1] <- -1

    d_chr <- 2 - 2 * sim
    d_chr[d_chr < 0] <- 0
    d_chr <- sqrt(d_chr)
    d_chr <- as.matrix(d_chr)
    if (!is.matrix(d_chr) || nrow(d_chr) != n_cells || ncol(d_chr) != n_cells) {
      stop("paper_innerproduct_mds produced invalid chromosome distance shape for ", chr, ": ",
           paste(dim(d_chr), collapse = "x"),
           " (expected ", n_cells, "x", n_cells, ")")
    }
    diag(d_chr) <- 0
    d_chr
  })

  d_blocks <- d_blocks[!vapply(d_blocks, is.null, logical(1))]
  if (length(d_blocks) == 0) {
    stop("paper_innerproduct_mds failed: no valid chromosome distance blocks")
  }

  arr <- simplify2array(d_blocks)
  if (length(dim(arr)) == 2L) {
    d_final <- arr
  } else {
    d_final <- apply(arr, c(1, 2), stats::median, na.rm = TRUE)
  }
  d_final <- as.matrix(d_final)
  if (!is.matrix(d_final) || nrow(d_final) != n_cells || ncol(d_final) != n_cells) {
    stop("paper_innerproduct_mds produced invalid distance shape: ",
         paste(dim(d_final), collapse = "x"),
         " (expected ", n_cells, "x", n_cells, ")")
  }

  d_final[!is.finite(d_final)] <- 0
  d_final <- (d_final + t(d_final)) / 2
  diag(d_final) <- 0
  d_final
}

embed_and_score <- function(Y_feature_by_cell, labels, method_name, pipeline, cfg, coverage_vec = NULL) {
  ari_uncorrected <- NA_real_

  if (identical(pipeline, "two_stage_pca_fixed")) {
    maxf <- cfg$two_stage_max_features_per_chr
    if (length(maxf) == 0 || is.null(maxf)) maxf <- 8000

    ts <- run_two_stage_pca(
      Y_feature_by_cell = Y_feature_by_cell,
      feature_names = rownames(Y_feature_by_cell),
      npc = cfg$npc,
      pcs_per_chr = cfg$schi_pcs_per_chr,
      max_features_per_chr = maxf
    )

    pcs0 <- ts$pcs_for_cluster
    k <- length(unique(labels))
    km0 <- kmeans(pcs0, centers = k, nstart = cfg$k_nstart)
    ari_uncorrected <- mclust::adjustedRandIndex(km0$cluster, as.integer(factor(labels)))

    pcs_cluster <- robust_clip_pcs(pcs0, clip_z = cfg$two_stage_pc_clip)
    if (ncol(pcs0) < 2) {
      emb <- cbind(pcs0[, 1], rep(0, nrow(pcs0)))
    } else {
      emb <- pcs0[, 1:2, drop = FALSE]
    }
    cluster_space <- pcs_cluster
  } else {
    X <- NULL
    get_X <- function() {
      if (is.null(X)) {
        X <<- prep_cell_feature_matrix(Y_feature_by_cell)
        if (ncol(X) < 2) stop("Too few non-constant features for embedding.")
      }
      X
    }

    if (identical(pipeline, "1d_pca") || identical(pipeline, "one_d_pca_visibility")) {
      max_bins <- cfg$one_d_max_bins
      if (length(max_bins) == 0 || is.null(max_bins) || !is.finite(max_bins) || max_bins < 10) {
        max_bins <- 5000L
      } else {
        max_bins <- as.integer(max_bins[1])
      }

      X1d <- build_cis_visibility_matrix(
        Y_feature_by_cell = Y_feature_by_cell,
        feature_names = rownames(Y_feature_by_cell),
        max_bins = max_bins
      )
      if (ncol(X1d) < 2) stop("1d_pca: too few non-constant bins")

      pca <- prcomp(X1d, center = TRUE, scale. = TRUE)
      npc_i <- min(cfg$npc, ncol(pca$x))
      pcs <- pca$x[, seq_len(npc_i), drop = FALSE]
      if (ncol(pcs) < 2) {
        emb <- cbind(pcs[, 1], rep(0, nrow(pcs)))
      } else {
        emb <- pcs[, 1:2, drop = FALSE]
      }
      cluster_space <- pcs
    } else
    if (identical(pipeline, "paper_innerproduct_mds")) {
      dmat <- build_innerproduct_distance_paper(
        Y_feature_by_cell = Y_feature_by_cell,
        feature_names = rownames(Y_feature_by_cell),
        n_strata = cfg$innerproduct_n_strata
      )
      dist_mat <- stats::as.dist(dmat)
      emb <- stats::cmdscale(dist_mat, k = 2)
      cluster_space <- emb
    } else
    if (identical(pipeline, "innerproduct_mds")) {
      Ymds <- Y_feature_by_cell
      if (!inherits(Ymds, "dgCMatrix")) {
        Ymds <- methods::as(Ymds, "dgCMatrix")
      }
      maxf_mds <- cfg$innerproduct_max_features
      if (length(maxf_mds) == 0 || is.null(maxf_mds) || !is.finite(maxf_mds) || maxf_mds < 2) {
        maxf_mds <- 10000L
      } else {
        maxf_mds <- as.integer(maxf_mds[1])
      }
      if (nrow(Ymds) > maxf_mds) {
        s <- Matrix::rowSums(Ymds)
        keep <- order(s, decreasing = TRUE)[seq_len(maxf_mds)]
        Ymds <- Ymds[keep, , drop = FALSE]
      }

      Xs <- Matrix::t(Ymds)
      rs <- Matrix::rowSums(Xs)
      rs[rs <= 0] <- 1
      Xs <- Matrix::Diagonal(x = 1 / rs) %*% Xs
      if (length(Xs@x) > 0) {
        Xs@x <- log1p(Xs@x * 1e6)
      }

      n_cells <- nrow(Xs)
      if (is.null(n_cells) || !is.finite(n_cells) || n_cells < 2) {
        stop("innerproduct_mds requires at least 2 cells; got ", n_cells)
      }

      sim_raw <- Matrix::tcrossprod(Xs)
      sim_vec <- as.numeric(sim_raw)
      if (length(sim_vec) != n_cells * n_cells) {
        stop("innerproduct_mds similarity length mismatch: ", length(sim_vec),
             " vs expected ", n_cells * n_cells)
      }
      sim <- matrix(sim_vec, nrow = n_cells, ncol = n_cells)

      norms <- sqrt(pmax(Matrix::rowSums(Xs^2), .Machine$double.eps))
      sim <- sim / outer(norms, norms)
      sim[!is.finite(sim)] <- 0
      sim <- pmin(1, pmax(-1, sim))

      dmat <- 1 - sim
      dmat[!is.finite(dmat)] <- 0
      dmat[dmat < 0] <- 0
      dmat <- base::as.matrix(dmat)
      if (!is.numeric(dmat)) storage.mode(dmat) <- "double"
      dd <- dim(dmat)
      if (is.null(dd) || length(dd) != 2L || dd[1] != n_cells || dd[2] != n_cells) {
        if (length(dmat) == n_cells * n_cells) {
          dmat <- matrix(as.numeric(dmat), nrow = n_cells, ncol = n_cells)
        }
      }
      dd <- dim(dmat)
      if (length(dd) != 2L || dd[1] != dd[2]) stop("innerproduct_mds produced non-square distance matrix: ", dd[1], "x", dd[2])
      diag(dmat) <- 0
      dist_mat <- stats::as.dist(dmat)
      emb <- stats::cmdscale(dist_mat, k = 2)
      cluster_space <- emb
    } else if (identical(pipeline, "pca_umap")) {
      X <- get_X()
      pca <- prcomp(X, center = TRUE, scale. = TRUE)
      npc_i <- min(cfg$npc, ncol(pca$x))
      pcs <- pca$x[, seq_len(npc_i), drop = FALSE]
      emb <- uwot::umap(
        pcs,
        n_neighbors = cfg$umap_n_neighbors,
        min_dist = cfg$umap_min_dist,
        metric = cfg$umap_metric,
        n_components = 2,
        verbose = FALSE,
        ret_model = FALSE
      )
      cluster_space <- pcs
    } else if (identical(pipeline, "schicluster_like")) {
      emb <- schicluster_like_embed(Y_feature_by_cell, cfg, force_pipeline = "pca_umap")
      if (is.null(emb)) {
        X <- get_X()
        pca <- prcomp(X, center = TRUE, scale. = TRUE)
        npc_i <- min(cfg$npc, ncol(pca$x))
        pcs <- pca$x[, seq_len(npc_i), drop = FALSE]
        emb <- uwot::umap(
          pcs,
          n_neighbors = cfg$umap_n_neighbors,
          min_dist = cfg$umap_min_dist,
          metric = cfg$umap_metric,
          n_components = 2,
          verbose = FALSE,
          ret_model = FALSE
        )
        cluster_space <- pcs
    } else {
        cluster_space <- emb
      }
    } else {
      stop("Unsupported pipeline: ", pipeline)
    }
  }

  k <- length(unique(labels))
  km <- kmeans(cluster_space, centers = k, nstart = cfg$k_nstart)
  ari <- mclust::adjustedRandIndex(km$cluster, as.integer(factor(labels)))

  tibble::tibble(
    X1 = emb[, 1],
    X2 = emb[, 2],
    cluster = factor(km$cluster),
    label = as.character(labels),
    method = method_name,
    pipeline = pipeline,
    ari = ari,
    ari_uncorrected = ari_uncorrected
  )
}













