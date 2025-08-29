#' Estimate a design-free AR(p) whitener from a run
#'
#' @description
#' Estimates an autoregressive whitener for a single `NeuroVec` run using the
#' global mean (optionally restricted to WM) and Yule–Walker AR fitting.
#' Returns the inverse Cholesky `Q`, the covariance `Sigma`, and AR
#' coefficients `phi`.
#'
#' @details
#' The time series is formed as the mean across voxels, optionally within a
#' provided white-matter mask. It is standardized, and AR(p) coefficients are
#' estimated via `stats::ar.yw`. The implied Toeplitz covariance `Sigma` is
#' regularized to positive definiteness (via a ridge and `Matrix::nearPD` if
#' needed). The whitener `Q = chol(Sigma)^{-1}` is returned for use in the
#' row metric. If `p = 0`, an identity whitener is returned.
#'
#' @param nv_run A `neuroim2::NeuroVec` (single run) with dimensions V x T.
#' @param wm_mask Optional `NeuroVol` to restrict the mean to WM voxels.
#' @param p Integer AR order (default 1). Use 0 to skip AR modeling.
#'
#' @return A list with components `Q` (inverse Cholesky), `Sigma` (covariance),
#'   and `phi` (AR coefficients of length `p`).
#'
#' @family temporal metrics
#' @seealso [estimate_ar_whitener_parcel()] for parcel-level AR whitening
#'
#' @examples
#' \donttest{
#' nv_run <- neuroim2::simulate_fmri(mask, n_time = 50, seed = 1)
#' arw   <- estimate_ar_whitener(nv_run, p = 1L)
#' }
#'
#' @export
estimate_ar_whitener <- function(nv_run, wm_mask = NULL, p = 1L) {
  stopifnot(inherits(nv_run, "NeuroVec"))
  X <- as.matrix(nv_run)     # V x T
  if (inherits(wm_mask, "NeuroVol")) {
    wv  <- neuroim2::values(wm_mask)
    widx <- which(wv > 0)
    y <- if (length(widx) > 0L) colMeans(X[widx, , drop = FALSE]) else colMeans(X)
  } else {
    y <- colMeans(X)
  }
  y <- as.numeric(base::scale(y, center = TRUE, scale = TRUE))
  Tlen <- length(y)

  # Handle p = 0 case (no AR modeling)
  if (p == 0) {
    Sigma <- Matrix::Diagonal(Tlen)
    Q <- Matrix::Diagonal(Tlen)
    phi <- numeric(0)
    return(list(Q = Q, Sigma = Sigma, phi = phi))
  }

  arf <- tryCatch(stats::ar.yw(y, order.max = p, aic = FALSE),
                  error = function(e) NULL)
  phi <- if (!is.null(arf) && length(arf$ar) >= 1L) as.numeric(arf$ar) else rep(0, p)

  acf <- stats::ARMAacf(ar = phi, lag.max = Tlen - 1L)
  Sigma <- toeplitz(as.numeric(acf))
  Sigma <- .make_pd(Sigma, eps = 1e-6)
  R     <- chol(Sigma)
  Q     <- solve(R)
  list(Q = Q, Sigma = Sigma, phi = phi)
}

.ar_corr <- function(phi, Tlen) {
  acf <- stats::ARMAacf(ar = phi, lag.max = Tlen - 1L)
  toeplitz(as.numeric(acf))
}

.make_pd <- function(S, eps = 1e-6) {
  S <- Matrix::forceSymmetric(S)
  ok <- TRUE
  tryCatch({ chol(S) }, error = function(e) ok <<- FALSE)
  if (!ok) {
    S <- S + Matrix::Diagonal(nrow(S)) * eps
    ok <- TRUE
    tryCatch({ chol(S) }, error = function(e) ok <<- FALSE)
    if (!ok && requireNamespace("Matrix", quietly = TRUE)) {
      S <- Matrix::nearPD(S, corr = TRUE)$mat
    }
  }
  Matrix::forceSymmetric(S)
}

#' Second-difference temporal smoothing penalty
#'
#' @description
#' Builds a T x T penalty matrix `I + lambda_t * D2' D2` where `D2` is the
#' second-difference operator. Promotes smoothness in time while keeping the
#' metric positive definite.
#'
#' @param Tlen Integer number of time points.
#' @param lambda_t Nonnegative weight on the roughness penalty (default 0.3).
#'
#' @return Symmetric positive definite `Matrix` of size `Tlen x Tlen`.
#'
#' @family temporal metrics
#' 
#' @examples
#' \donttest{
#' H <- make_temporal_penalty(30, lambda_t = 0.5)
#' }
#'
#' @export
make_temporal_penalty <- function(Tlen, lambda_t = 0.3) {
  # Create second difference matrix explicitly
  # D2 is (T-2) x T matrix representing second differences
  if (Tlen <= 2) {
    # For very short time series, just return identity
    return(Matrix::Diagonal(Tlen))
  }
  
  # Build second difference matrix manually
  i <- rep(seq_len(Tlen - 2), each = 3)
  j <- as.vector(sapply(seq_len(Tlen - 2), function(k) k:(k+2)))
  x <- rep(c(1, -2, 1), Tlen - 2)
  D2 <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(Tlen - 2, Tlen))
  
  Matrix::Diagonal(Tlen) + lambda_t * Matrix::forceSymmetric(Matrix::t(D2) %*% D2)
}

#' Frame-wise robust weights from FD and DVARS
#'
#' @description
#' Converts motion summaries (FD) and DVARS into diagonal frame weights in (0,1],
#' down-weighting high-motion or high-DVARS frames. If both inputs are `NULL`,
#' returns a `1x1` identity (caller may expand to `T x T`).
#'
#' @details
#' ## Framewise Displacement (FD)
#' 
#' FD quantifies head motion between consecutive fMRI volumes as defined by
#' Power et al. (2012). It is calculated as the sum of absolute values of the
#' temporal derivatives of the six rigid-body motion parameters (3 translations
#' and 3 rotations):
#' 
#' `FD = |Δx| + |Δy| + |Δz| + |Δα| + |Δβ| + |Δγ|`
#' 
#' where rotational parameters (α, β, γ) are converted to millimeters by 
#' calculating arc length displacement on a sphere with radius 50mm. This
#' conversion ensures all parameters are on the same scale. FD values represent
#' the total displacement in mm from one volume to the next.
#' 
#' ## DVARS (Derivative of VARS)
#' 
#' DVARS measures the rate of change of BOLD signal intensity across the whole
#' brain between consecutive time points (Smyser et al., 2010; Power et al., 2012).
#' It is calculated as the root mean square of the temporal derivative:
#' 
#' `DVARS_t = sqrt(mean((x_t - x_{t-1})^2))`
#' 
#' where x represents voxel intensities. DVARS is sensitive to both motion
#' artifacts and physiological noise. High DVARS values indicate sudden global
#' signal changes that often correlate with head motion.
#' 
#' ## Weight Calculation
#' 
#' This function converts these metrics to frame weights:
#' - FD weights: `w_FD = 1 / (1 + max(0, FD - fd_thresh))`
#' - DVARS weights: `w_DVARS = 1 / (1 + max(0, z - dvars_z))` where z is the z-score
#' - Final weight: `w = max(0.001, w_FD * w_DVARS)`
#' 
#' The resulting diagonal weight matrix down-weights problematic frames while
#' preserving data continuity, as complete censoring can introduce discontinuities
#' in the time series.
#'
#' @param FD Optional numeric vector of framewise displacement values in mm.
#'   Typically computed from motion correction parameters during preprocessing.
#' @param DVARS Optional numeric vector of DVARS values (signal RMS derivative).
#'   Typically computed from preprocessed BOLD time series.
#' @param fd_thresh Threshold above which FD begins to down-weight (default 0.5mm).
#'   Common choices: 0.2mm (strict), 0.5mm (standard), 1.0mm (lenient).
#' @param dvars_z Z-score threshold for DVARS down-weighting (default 3.0).
#'   Frames with DVARS z-scores above this are down-weighted; 2.5-3.5 typical.
#'
#' @return A diagonal `Matrix` of weights.
#'
#' @family temporal metrics
#' 
#' @references
#' Power JD, Barnes KA, Snyder AZ, Schlaggar BL, Petersen SE (2012).
#' "Spurious but systematic correlations in functional connectivity MRI
#' networks arise from subject motion." NeuroImage, 59(3), 2142-2154.
#' 
#' Smyser CD, Inder TE, Shimony JS, Hill JE, Degnan AJ, Snyder AZ, Neil JJ (2010).
#' "Longitudinal analysis of neural network development in preterm infants."
#' Cerebral Cortex, 20(12), 2852-2862.
#' 
#' @examples
#' # Create frame weights from motion parameters
#' W <- make_frame_weights(FD = runif(30, 0, 0.5), DVARS = runif(30, 0.5, 1.5))
#' diag(W)[1:5]  # First 5 frame weights
#'
#' @export
make_frame_weights <- function(FD = NULL, DVARS = NULL,
                               fd_thresh = 0.5, dvars_z = 3.0) {
  # Handle case when both are NULL - return 1x1 identity
  if (is.null(FD) && is.null(DVARS)) {
    return(Matrix::Diagonal(1))
  }
  
  if (is.null(FD))    FD    <- rep(0, length(DVARS))
  if (is.null(DVARS)) DVARS <- rep(0, length(FD))
  L <- max(length(FD), length(DVARS))
  FD <- as.numeric(FD)[seq_len(L)]
  DVARS <- as.numeric(DVARS)[seq_len(L)]
  fdw  <- 1 / (1 + pmax(0, FD - fd_thresh))
  z    <- base::scale(DVARS)
  dvrw <- 1 / (1 + pmax(0, as.numeric(z) - dvars_z))
  w    <- pmax(1e-3, fdw * dvrw)
  Matrix::Diagonal(x = w)
}

#' Build temporal covariance metric for generalized PCA
#'
#' @description
#' Constructs a temporal covariance metric (row metric M) that encodes temporal
#' dependencies, smoothness constraints, and data quality weights for a single
#' fMRI run. This metric is used as the row-wise (temporal) metric in generalized
#' PCA to properly account for temporal structure in fMRI data.
#'
#' @details
#' ## Purpose in Generalized PCA
#' 
#' In standard PCA, all time points are treated as independent and equally
#' reliable. This is inappropriate for fMRI data where:
#' - Consecutive time points are correlated (temporal autocorrelation)
#' - Some time points are corrupted by motion artifacts
#' - Neural signals vary smoothly over time
#' 
#' This function builds a metric M that encodes these properties, allowing
#' generalized PCA to find components that respect the temporal structure
#' of fMRI data.
#' 
#' ## Mathematical Formulation
#' 
#' The temporal metric combines three key components:
#' 
#' `M = W^{1/2} Q' H Q W^{1/2} + ridge * I`
#' 
#' where:
#' - `Q`: AR(p) whitening matrix that removes temporal autocorrelation
#' - `H`: Temporal smoothness penalty (second-difference operator)
#' - `W`: Diagonal weights based on motion quality (FD/DVARS)
#' - `ridge`: Small regularization for numerical stability
#' 
#' ## Components Explained
#' 
#' **AR Whitening (Q)**: fMRI time series exhibit strong autocorrelation
#' (typically 0.3-0.5 at lag 1). The AR whitener transforms the data to make
#' time points approximately independent, preventing overrepresentation of
#' slow temporal drifts in the components.
#' 
#' **Temporal Smoothness (H)**: Neural signals change smoothly over time.
#' The penalty `H = I + λ_t * D2' D2` (where D2 is the second-difference
#' operator) encourages components with smooth temporal profiles while
#' allowing for dynamic changes.
#' 
#' **Motion Weighting (W)**: Time points with high motion (FD) or global
#' signal changes (DVARS) are down-weighted rather than discarded. This
#' soft censoring approach maintains temporal continuity while reducing
#' artifact influence. See `make_frame_weights()` for details.
#' 
#' ## Usage Scenarios
#' 
#' - **Standard**: `p=1, lambda_t=0.3` with both FD and DVARS
#' - **High motion data**: Increase lambda_t (0.5-0.7) for more smoothing
#' - **Task fMRI**: Consider p=0 if modeling task regressors separately
#' - **Resting state**: Always use p≥1 to account for autocorrelation
#'
#' @param nv_run A `neuroim2::NeuroVec` (single run) containing fMRI time series.
#' @param wm_mask Optional `NeuroVol` white matter mask. When provided, AR
#'   parameters are estimated from the mean white matter signal, which typically
#'   has less neural signal and better represents noise autocorrelation.
#' @param p Integer AR order (default 1). Common choices:
#'   - 0: No AR modeling (independent time points)
#'   - 1: First-order autoregression (standard, captures ~70% of autocorrelation)
#'   - 2-3: Higher orders for data with stronger temporal dependencies
#' @param FD Optional numeric vector of framewise displacement values (mm).
#'   Length must equal the number of time points. Typically computed during
#'   motion correction preprocessing.
#' @param DVARS Optional numeric vector of DVARS values (signal derivative RMS).
#'   Length must equal the number of time points. Typically computed from
#'   preprocessed BOLD signal.
#' @param lambda_t Temporal smoothness penalty weight (default 0.3). Range 0-1:
#'   - 0: No temporal smoothing
#'   - 0.3: Standard smoothing (default)
#'   - 0.5-0.7: Strong smoothing for high-motion or noisy data
#'   - 1.0: Maximum smoothing
#' @param ridge Small positive value added to diagonal for numerical stability
#'   (default 1e-6). Increase if encountering singular matrix errors.
#'
#' @return A symmetric positive definite `Matrix` of size `T x T` encoding
#'   temporal covariance structure for use in generalized PCA.
#'
#' @family temporal metrics
#' @seealso 
#' - [build_temporal_metric_parcel()] for parcel-level temporal metrics
#' - [estimate_ar_whitener()] for AR parameter estimation details
#' - [make_temporal_penalty()] for smoothness penalty construction
#' - [make_frame_weights()] for motion-based weighting
#' 
#' @examples
#' \donttest{
#' # Basic usage with AR(1) whitening
#' M <- build_temporal_metric(nv_run, p = 1L)
#' 
#' # With motion correction
#' M <- build_temporal_metric(nv_run, wm_mask = wm_vol, 
#'                            FD = fd_values, DVARS = dvars_values,
#'                            lambda_t = 0.5)
#' 
#' # No AR modeling for task fMRI with explicit task regressors
#' M <- build_temporal_metric(nv_run, p = 0, lambda_t = 0.3)
#' }
#'
#' @export
build_temporal_metric <- function(nv_run, wm_mask = NULL, p = 1L,
                                  FD = NULL, DVARS = NULL,
                                  lambda_t = 0.3, ridge = 1e-6) {
  Tlen <- dim(nv_run)[4]
  Q    <- estimate_ar_whitener(nv_run, wm_mask, p = p)$Q
  H    <- make_temporal_penalty(Tlen, lambda_t)
  W    <- make_frame_weights(FD, DVARS)
  
  # If no motion parameters provided, make W match Tlen dimensions
  if (nrow(W) == 1 && ncol(W) == 1) {
    W <- Matrix::Diagonal(Tlen)
  }
  
  M    <- Matrix::forceSymmetric((W^(1/2)) %*% (Matrix::t(Q) %*% H %*% Q) %*% (W^(1/2)))
  M + Matrix::Diagonal(Tlen) * ridge
}
