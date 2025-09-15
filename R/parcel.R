#' Parcel centroids from a label volume
#'
#' @description
#' Computes centroid coordinates (i,j,k) for labeled parcels in a
#' `NeuroVol` label image. Labels with value 0 are ignored.
#'
#' @details
#' For each parcel label, voxel coordinates are averaged to obtain a centroid.
#' Optionally restrict to a provided vector of `labels`.
#'
#' @param parc_vol `NeuroVol` of parcel labels (integer; 0 treated as background).
#' @param labels Optional integer vector of labels to include.
#'
#' @return A list with `coords` (P x 3 matrix) and `labels` (length P).
#'
#' @family parcel utilities
#'
#' @examples
#' \donttest{
#' pc <- parcel_centroids_from_labels(parc_vol)
#' }
#'
#' @export
parcel_centroids_from_labels <- function(parc_vol, labels = NULL) {
  stopifnot(inherits(parc_vol, "NeuroVol"))
  lab_all <- as.integer(neuroim2::values(parc_vol))
  labs <- unique(lab_all[lab_all > 0L]); labs <- sort(labs)
  if (length(labs) == 0) {
    stop("No non-zero labels found in parcellation volume")
  }
  if (!is.null(labels)) labs <- as.integer(labels)

  dims <- dim(parc_vol)
  Vx <- dims[1]; Vy <- dims[2]; Vz <- if (length(dims) >= 3) dims[3] else 1L

  idx <- which(lab_all > 0L)
  lin <- idx - 1L
  ii <- (lin %% Vx) + 1L
  jj <- ((lin %/% Vx) %% Vy) + 1L
  kk <- (lin %/% (Vx * Vy)) + 1L
  coords_all <- cbind(ii, jj, kk)
  labs_all   <- lab_all[idx]

  P <- length(labs)
  C <- matrix(NA_real_, nrow = P, ncol = 3L)
  for (p in seq_len(P)) {
    L <- labs[p]
    sel <- labs_all == L
    if (!any(sel)) stop("Label ", L, " not found in parc_vol.")
    C[p, ] <- colMeans(coords_all[sel, , drop = FALSE])
  }
  colnames(C) <- c("i","j","k")
  list(coords = C, labels = labs)
}

#' Aggregate GM/WM/CSF tissue probabilities to parcels
#'
#' @description
#' Computes mean tissue probabilities within each parcel label to obtain
#' parcel-level summaries suitable for parcel metrics.
#'
#' @param gm_vol,wm_vol,csf_vol `NeuroVol` tissue maps.
#' @param parc_vol `NeuroVol` of parcel labels.
#' @param labels Optional integer vector of labels to include/order.
#'
#' @return A list with numeric vectors `gm`, `wm`, `csf`, and `labels`.
#'
#' @family parcel utilities
#'
#' @examples
#' \donttest{
#' # Example requires tissue maps and parcellation
#' # See vignette for complete example
#' }
#'
#' @export
aggregate_tissue_to_parcels <- function(gm_vol, wm_vol, csf_vol, parc_vol, labels = NULL) {
  stopifnot(inherits(gm_vol, "NeuroVol"), inherits(wm_vol, "NeuroVol"),
            inherits(csf_vol, "NeuroVol"), inherits(parc_vol, "NeuroVol"))
  sp <- neuroim2::space(parc_vol)
  if (!identical(neuroim2::space(gm_vol), sp) ||
      !identical(neuroim2::space(wm_vol), sp) ||
      !identical(neuroim2::space(csf_vol), sp)) {
    stop("Tissue maps and parc_vol must share the same NeuroSpace.")
  }
  lab_all <- as.integer(neuroim2::values(parc_vol))
  labs <- unique(lab_all[lab_all > 0L]); labs <- sort(labs)
  if (!is.null(labels)) labs <- as.integer(labels)

  idx <- which(lab_all > 0L)
  labs_idx <- lab_all[idx]
  gm  <- as.numeric(neuroim2::values(gm_vol))[idx]
  wm  <- as.numeric(neuroim2::values(wm_vol))[idx]
  csf <- as.numeric(neuroim2::values(csf_vol))[idx]

  gm_p  <- tapply(gm,  labs_idx, function(x) mean(x, na.rm = TRUE))
  wm_p  <- tapply(wm,  labs_idx, function(x) mean(x, na.rm = TRUE))
  csf_p <- tapply(csf, labs_idx, function(x) mean(x, na.rm = TRUE))

  gm_p  <- as.numeric(gm_p[as.character(labs)])
  wm_p  <- as.numeric(wm_p[as.character(labs)])
  csf_p <- as.numeric(csf_p[as.character(labs)])

  list(gm = gm_p, wm = wm_p, csf = csf_p, labels = labs)
}

#' Build a parcel-level Laplacian
#'
#' @description
#' Constructs a parcel adjacency-based Laplacian either from a supplied adjacency
#' matrix, parcel centroids (k–NN), or a label volume. Normalization and spectral
#' scaling mirror the voxel-level Laplacian.
#'
#' @param adjacency Optional symmetric adjacency/weight matrix (P x P).
#' @param parcel_coords Optional P x 3 matrix of parcel centroids.
#' @param parc_vol Optional `NeuroVol` of parcel labels (used to compute centroids).
#' @param k,sigma KNN neighbors and Gaussian bandwidth.
#' @param normalized Logical; return normalized Laplacian.
#'
#' @return Sparse symmetric `Matrix` (P x P) Laplacian.
#'
#' @family spatial metrics
#' @seealso [make_laplacian()] for voxel-level Laplacian construction
#'
#' @examples
#' \donttest{
#' Lp <- make_parcel_laplacian(parc_vol = parc, k = 6, sigma = 2)
#' }
#'
#' @export
make_parcel_laplacian <- function(adjacency = NULL, parcel_coords = NULL, parc_vol = NULL,
                                  k = 6, sigma = 2.5, normalized = TRUE) {
  if (!is.null(adjacency)) {
    W <- Matrix::forceSymmetric(adjacency)
  } else {
    if (is.null(parcel_coords)) {
      if (is.null(parc_vol)) stop("Provide either adjacency, parcel_coords, or parc_vol.")
      pc <- parcel_centroids_from_labels(parc_vol)
      parcel_coords <- pc$coords
    }
    W <- .knn_gaussian(parcel_coords, k = k, sigma = sigma, symmetric = TRUE)
  }

  d <- Matrix::rowSums(W)
  D <- Matrix::Diagonal(x = as.numeric(pmax(d, 1e-12)))
  if (normalized) {
    Dm <- Matrix::Diagonal(x = 1 / sqrt(D@x))
    I  <- Matrix::Diagonal(nrow(W))
    L  <- Matrix::forceSymmetric(I - Dm %*% W %*% Dm)
  } else {
    L  <- Matrix::forceSymmetric(D - W)
  }

  lam_max <- if (requireNamespace("RSpectra", quietly = TRUE) && nrow(L) >= 100) {
    as.numeric(RSpectra::eigs_sym(L, 1, which = "LM")$values)
  } else {
    max(Re(eigen(as.matrix(L), only.values = TRUE)$values))
  }
  if (is.finite(lam_max) && lam_max > 0) L <- L / lam_max
  L
}

#' Parcel-level tSNR across runs (ClusteredNeuroVec)
#'
#' @description
#' Computes parcel-wise temporal signal-to-noise ratio across `ClusteredNeuroVec`
#' runs by aggregating means and variances over time.
#'
#' @details
#' For each parcel, the mean and variance are computed over concatenated time
#' points across all runs. tSNR is calculated as `mean / sd`, with both
#' numerator and denominator lower-bounded at `1e-6` to avoid division by zero
#' or degenerate weights. This is the parcel-level analogue of `compute_tsnr()`
#' for voxel data.
#'
#' @param cnv_list List of `neuroim2::ClusteredNeuroVec` runs.
#'
#' @return Numeric vector of length P with parcel tSNR values.
#'
#' @family spatial metrics
#' @seealso [compute_tsnr()] for voxel-level tSNR computation
#'
#' @examples
#' \donttest{
#' # Compute parcel tSNR across multiple runs
#' tsnr_p <- compute_tsnr_parcel(cnv_list)
#' }
#'
#' @export
compute_tsnr_parcel <- function(cnv_list) {
  stopifnot(length(cnv_list) >= 1L)
  X0 <- as.matrix(cnv_list[[1L]])
  P  <- nrow(X0)
  s1 <- numeric(P); s2 <- numeric(P); Ttot <- 0L
  for (cnv in cnv_list) {
    X <- as.matrix(cnv)   # P x T
    s1 <- s1 + rowSums(X)
    s2 <- s2 + rowSums(X^2)
    Ttot <- Ttot + ncol(X)
  }
  mu  <- s1 / Ttot
  var <- pmax(0, s2 / Ttot - mu^2)
  sdv <- sqrt(var)
  pmax(1e-6, mu) / pmax(1e-6, sdv)
}

#' Estimate AR(p) whitener from parcel-level data
#'
#' @description
#' Parcel analogue of `estimate_ar_whitener()` that uses the parcel mean or
#' optional `wm_parcels` subset to estimate AR(p) via Yule–Walker.
#'
#' @details
#' The time series is formed as the mean across parcels, optionally restricted
#' to white-matter parcels if `wm_parcels` indices are provided. The series is
#' standardized and AR(p) coefficients are estimated via `stats::ar.yw`. The
#' implied Toeplitz covariance matrix is regularized to ensure positive
#' definiteness. Returns the whitener `Q = chol(Sigma)^{-1}` for use in the
#' row metric construction.
#'
#' @param cnv_run A `ClusteredNeuroVec` run.
#' @param wm_parcels Optional integer indices of WM parcels for more robust
#'   AR estimation.
#' @param p Integer AR order (default 1). Use higher orders for data with
#'   stronger temporal autocorrelation.
#'
#' @return List with components:
#'   \describe{
#'     \item{Q}{Inverse Cholesky factor (whitening matrix)}
#'     \item{Sigma}{Regularized AR covariance matrix}
#'     \item{phi}{AR coefficients of length p}
#'   }
#'
#' @family temporal metrics
#' @seealso [estimate_ar_whitener()] for voxel-level AR whitening
#'
#' @examples
#' \donttest{
#' # Estimate AR(1) whitener using WM parcels
#' arw <- estimate_ar_whitener_parcel(cnv_run, wm_parcels = c(1,2,3), p = 1)
#' }
#'
#' @export
estimate_ar_whitener_parcel <- function(cnv_run, wm_parcels = NULL, p = 1L) {
  stopifnot(inherits(cnv_run, "ClusteredNeuroVec"))
  X <- as.matrix(cnv_run)  # P x T
  y <- if (!is.null(wm_parcels) && length(wm_parcels) > 0L) colMeans(X[wm_parcels, , drop = FALSE]) else colMeans(X)
  y <- as.numeric(scale(y, center = TRUE, scale = TRUE))
  Tlen <- length(y)

  arf <- tryCatch(stats::ar.yw(y, order.max = p, aic = FALSE), error = function(e) NULL)
  phi <- if (!is.null(arf) && length(arf$ar) >= 1L) as.numeric(arf$ar) else rep(0, p)

  Sigma <- stats::toeplitz(stats::ARMAacf(ar = phi, lag.max = Tlen - 1L))
  Sigma <- .make_pd(Sigma, eps = 1e-6)
  R     <- chol(Sigma); Q <- solve(R)
  list(Q = Q, Sigma = Sigma, phi = phi)
}

#' Build temporal covariance metric for parcellated fMRI data
#'
#' @description
#' Constructs a temporal covariance metric (row metric M) for parcellated fMRI
#' data in `ClusteredNeuroVec` format. This is the parcel-level analogue of
#' `build_temporal_metric()` for voxel data.
#'
#' @details
#' The metric is `M = W^{1/2} Q' H Q W^{1/2} + ridge * I`, where:
#' \itemize{
#'   \item `Q` whitens the temporal covariance from AR(p) estimation
#'   \item `H` is the second-difference temporal smoothness penalty
#'   \item `W` is diagonal frame weighting from motion parameters
#' }
#' This parcel-level implementation uses parcel-averaged time series for AR
#' estimation, providing computational efficiency for high-dimensional data.
#'
#' @param cnv_run A `ClusteredNeuroVec` run.
#' @param wm_parcels Optional integer indices of WM parcels for AR estimation.
#' @param p AR order (default 1). Higher orders model longer-range dependencies.
#' @param FD,DVARS Optional motion metrics for the run. FD is framewise 
#'   displacement (mm), DVARS is temporal derivative of signal variance.
#'   See \code{\link{make_frame_weights}} for detailed definitions.
#' @param lambda_t Temporal penalty weight (default 0.3). Higher values enforce
#'   more temporal smoothness.
#' @param ridge Small ridge added to ensure positive definiteness (default 1e-6).
#'
#' @return Symmetric positive definite `Matrix` (T x T).
#'
#' @family temporal metrics
#' @seealso [build_temporal_metric()] for voxel-level temporal metric construction
#'
#' @examples
#' \donttest{
#' # Build temporal metric with AR(1) whitening and motion regressors
#' M <- build_temporal_metric_parcel(cnv_run, wm_parcels = c(1,2,3),
#'                                   FD = fd_vals, DVARS = dvars_vals)
#' }
#'
#' @export
build_temporal_metric_parcel <- function(cnv_run, wm_parcels = NULL, p = 1L,
                                         FD = NULL, DVARS = NULL,
                                         lambda_t = 0.3, ridge = 1e-6) {
  Tlen <- dim(cnv_run)[2]
  Q    <- estimate_ar_whitener_parcel(cnv_run, wm_parcels = wm_parcels, p = p)$Q
  H    <- make_temporal_penalty(Tlen, lambda_t)
  W    <- make_frame_weights(FD, DVARS)
  M    <- Matrix::forceSymmetric((W^(1/2)) %*% (Matrix::t(Q) %*% H %*% Q) %*% (W^(1/2)))
  M + Matrix::Diagonal(Tlen) * ridge
}

#' Build parcel-level spatial column metric A from aggregated tissues and Laplacian
#'
#' @description
#' Forms a parcel-level spatial metric `A` analogous to voxel-level `build_spatial_metric`,
#' using aggregated tissue probabilities and optional parcel tSNR.
#'
#' @details
#' Constructs weights `w = gm^alpha * tsnr^beta * (wm + csf)^{-gamma}` for each
#' parcel, then forms `A = diag(sqrt(w)) * (I + lambda_s * Lp) * diag(sqrt(w)) + tau * I`.
#' This encodes both spatial smoothness (via the Laplacian) and parcel-wise
#' quality weights derived from tissue composition and temporal SNR. The result
#' is a positive definite metric suitable for generalized PCA.
#'
#' @param gm_p,wm_p,csf_p Numeric vectors of length P with parcel tissue values
#'   (typically means of voxel-wise tissue probabilities within each parcel).
#' @param Lp Parcel Laplacian (P x P) from `make_parcel_laplacian()`.
#' @param tsnr_p Optional parcel tSNR vector from `compute_tsnr_parcel()`.
#' @param alpha Exponent for GM weighting (default 1.0). Higher values give more
#'   weight to gray matter parcels.
#' @param beta Exponent for tSNR weighting (default 0.5). Controls influence of
#'   temporal signal quality.
#' @param gamma Exponent for (WM+CSF) down-weighting (default 1.0). Higher values
#'   penalize non-gray matter more strongly.
#' @param lambda_s Spatial regularization weight (default 0.5). Controls the
#'   strength of spatial smoothing via the Laplacian.
#' @param tau Small ridge on the diagonal (default 1e-6) for numerical stability.
#'
#' @return Symmetric positive definite `Matrix` (P x P).
#'
#' @family spatial metrics
#' @seealso [build_spatial_metric()] for voxel-level column metric construction
#'
#' @examples
#' \donttest{
#' # Build parcel column metric with tissue weights and spatial smoothing
#' gm_p <- c(0.8, 0.6, 0.2)
#' wm_p <- c(0.15, 0.25, 0.6)
#' csf_p <- c(0.05, 0.15, 0.2)
#' tsnr_vals <- c(50, 40, 30)
#' adj <- Matrix::Matrix(matrix(c(0, 1, 0,
#'                                1, 0, 1,
#'                                0, 1, 0), nrow = 3, byrow = TRUE))
#' deg <- Matrix::Diagonal(x = c(1, 2, 1))
#' Lp <- deg - adj
#' A <- build_spatial_metric_parcel(gm_p, wm_p, csf_p, Lp, tsnr_p = tsnr_vals,
#'                                  lambda_s = 0.5)
#' }
#'
#' @export
build_spatial_metric_parcel <- function(gm_p, wm_p, csf_p, Lp,
                                tsnr_p = NULL,
                                alpha = 1.0, beta = 0.5, gamma = 1.0,
                                lambda_s = 0.5, tau = 1e-6) {
  stopifnot(length(gm_p) == length(wm_p), length(wm_p) == length(csf_p))
  P <- length(gm_p)
  if (is.null(tsnr_p)) tsnr_p <- rep(1, P)
  w <- (pmax(gm_p, 1e-6)^alpha) * (pmax(tsnr_p, 1e-6)^beta) * ((pmax(wm_p + csf_p, 1e-6))^(-gamma))
  S <- Matrix::Diagonal(P) + lambda_s * Lp
  A <- Matrix::forceSymmetric((Matrix::Diagonal(x = sqrt(w)) %*% S %*% Matrix::Diagonal(x = sqrt(w))))
  A + Matrix::Diagonal(P) * tau
}

#' Fit generalized PCA at the parcel level (ClusteredNeuroVec)
#'
#' @description
#' Builds a parcel Laplacian and parcel-level `A`/`M` metrics then fits genpca
#' on stacked runs. See voxel-level `fit_subject_genpca()` for an overview.
#'
#' @details
#' This is the parcel-level analogue of `fit_subject_genpca()`, operating on
#' `ClusteredNeuroVec` data where voxels have been aggregated into parcels.
#' The workflow is:
#' \enumerate{
#'   \item Build parcel graph Laplacian from adjacency, coordinates, or labels
#'   \item Compute parcel tSNR across runs
#'   \item Aggregate tissue maps to parcels (if not provided)
#'   \item Construct column metric A with spatial and tissue weights
#'   \item Build per-run row metrics M with AR whitening and motion weights
#'   \item Stack runs temporally and fit generalized PCA
#' }
#' Parcellation reduces dimensionality and computation while preserving
#' spatial structure through the parcel graph.
#'
#' @param cnv_list List of `ClusteredNeuroVec` runs.
#' @param parc_vol,parcel_coords,adjacency Ways to define the parcel graph:
#'   provide exactly one.
#' @param gm_vol,wm_vol,csf_vol Optional tissue maps used if `gm_p/wm_p/csf_p`
#'   are not provided.
#' @param gm_p,wm_p,csf_p Optional parcel-level tissue summaries (avoids
#'   recomputation).
#' @param wm_parcels Optional WM parcel indices for AR whitener estimation.
#' @param FD_list,DVARS_list Optional lists of per-run motion metrics.
#'   Each list should contain numeric vectors (one per run in cnv_list).
#'   FD is framewise displacement (mm), DVARS is temporal derivative of
#'   signal variance. Used to down-weight high-motion frames. See
#'   \code{\link{make_frame_weights}} for detailed definitions of these metrics.
#' @param k Number of components to extract (default 20).
#' @param p_ar AR order for temporal whitening (default 1).
#' @param lambda_s Spatial regularization weight (default 0.5).
#' @param lambda_t Temporal regularization weight (default 0.3).
#' @param knn_k Number of nearest neighbors for parcel graph (default 6).
#' @param knn_sigma Gaussian kernel bandwidth for parcel graph (default 2.5).
#'
#' @return A `genpca` object with parcel-level loadings.
#'
#' @family main fitting
#' @seealso [fit_subject_genpca()] for voxel-level analysis
#'
#' @examples
#' \donttest{
#' # Fit genpca on parcellated data with tissue weights
#' fit <- fit_subject_genpca_parcel(cnv_list, parc_vol = parcellation,
#'                                  gm_vol = gm, wm_vol = wm, csf_vol = csf,
#'                                  k = 20)
#' }
#'
#' @export
fit_subject_genpca_parcel <- function(cnv_list,
                                      parc_vol = NULL, parcel_coords = NULL, adjacency = NULL,
                                      gm_vol = NULL, wm_vol = NULL, csf_vol = NULL,
                                      gm_p = NULL, wm_p = NULL, csf_p = NULL,
                                      wm_parcels = NULL,
                                      FD_list = NULL, DVARS_list = NULL,
                                      k = 20, p_ar = 1L,
                                      lambda_s = 0.5, lambda_t = 0.3,
                                      knn_k = 6, knn_sigma = 2.5) {
  stopifnot(length(cnv_list) >= 1L)
  X0 <- as.matrix(cnv_list[[1L]])
  P  <- nrow(X0)

  Lp <- if (!is.null(adjacency)) {
    make_parcel_laplacian(adjacency = adjacency, normalized = TRUE)
  } else if (!is.null(parcel_coords)) {
    make_parcel_laplacian(parcel_coords = parcel_coords, k = knn_k, sigma = knn_sigma, normalized = TRUE)
  } else if (!is.null(parc_vol)) {
    make_parcel_laplacian(parc_vol = parc_vol, k = knn_k, sigma = knn_sigma, normalized = TRUE)
  } else {
    stop("Provide one of: adjacency, parcel_coords, or parc_vol to build the parcel Laplacian.")
  }

  tsnr_p <- compute_tsnr_parcel(cnv_list)

  if (is.null(gm_p) || is.null(wm_p) || is.null(csf_p)) {
    if (is.null(parc_vol) || is.null(gm_vol) || is.null(wm_vol) || is.null(csf_vol)) {
      stop("To aggregate tissues to parcels, provide parc_vol and gm_vol/wm_vol/csf_vol, or pass gm_p/wm_p/csf_p.")
    }
    agg <- aggregate_tissue_to_parcels(gm_vol, wm_vol, csf_vol, parc_vol)
    gm_p <- agg$gm_p; wm_p <- agg$wm_p; csf_p <- agg$csf_p
    if (length(gm_p) != P) stop("Number of parcels in parc_vol/tissue maps does not match ClusteredNeuroVec (P).")
  }

  A <- build_spatial_metric_parcel(gm_p, wm_p, csf_p, Lp, tsnr_p = tsnr_p, lambda_s = lambda_s)

  X_rows <- list(); M_blks <- list()
  for (r in seq_along(cnv_list)) {
    cnv <- cnv_list[[r]]
    Xpr <- t(as.matrix(cnv))                   # T_r x P
    FD_r    <- if (!is.null(FD_list))    FD_list[[r]]    else NULL
    DVARS_r <- if (!is.null(DVARS_list)) DVARS_list[[r]] else NULL
    Mr <- build_temporal_metric_parcel(cnv, wm_parcels = wm_parcels, p = p_ar,
                                       FD = FD_r, DVARS = DVARS_r,
                                       lambda_t = lambda_t, ridge = 1e-6)
    X_rows[[r]] <- Xpr
    M_blks[[r]] <- Mr
  }
  X_stack <- do.call(rbind, X_rows)
  M_blk   <- Matrix::bdiag(M_blks)

  genpca::genpca(X_stack, A = A, M = M_blk, ncomp = k, preproc = multivarious::center())
}
