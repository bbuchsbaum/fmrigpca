#' Form a whitened matrix for genpca geometry
#'
#' @description
#' Applies square-root scalings of the row and column metrics to obtain the
#' whitened matrix `W = R_M^T X R_A`, where `A = R_A^T R_A` and `M = R_M^T R_M`.
#'
#' @param X Numeric matrix of size T x V (time x voxels/parcels).
#' @param A Column metric (V x V), symmetric positive definite.
#' @param M Row metric (T x T), symmetric positive definite.
#'
#' @return The whitened matrix `W` as a `Matrix`.
#'
#' @family subspace analysis
#' @seealso [whitened_svd()] for computing SVD of the whitened matrix
#' 
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' A <- Matrix::Diagonal(ncol(X))
#' M <- Matrix::Diagonal(nrow(X))
#' W <- whitened_matrix(X, A, M)
#' }
#'
#' @export
whitened_matrix <- function(X, A, M) {
  stopifnot(nrow(X) == nrow(M), ncol(X) == nrow(A))
  RA <- chol(A)  # A = RA' RA
  RM <- chol(M)  # M = RM' RM
  (Matrix::t(RM)) %*% X %*% (RA)
}

#' Truncated SVD of the whitened matrix
#'
#' @description
#' Computes singular values (optionally truncated) of the whitened matrix. Uses
#' `RSpectra` when available for efficiency in the truncated case.
#'
#' @param X,A,M As in `whitened_matrix()`.
#' @param k Optional integer number of singular values to compute.
#'
#' @return A list with `d` (singular values) and `u`,`v` set to `NULL`.
#'
#' @family subspace analysis
#' @seealso [whitened_matrix()] for the underlying whitening operation
#' 
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(300), nrow = 30, ncol = 10)
#' A <- Matrix::Diagonal(ncol(X))
#' M <- Matrix::Diagonal(nrow(X))
#' ws <- whitened_svd(X, A, M, k = 5)
#' }
#'
#' @export
whitened_svd <- function(X, A, M, k = NULL) {
  W <- whitened_matrix(X, A, M)
  if (is.null(k)) {
    sv <- svd(as.matrix(W), nu = 0, nv = 0)
    return(list(d = sv$d, u = NULL, v = NULL))
  }
  # Truncated SVD operating directly on W (avoids forming the dense V x V
  # crossproduct, which is infeasible for voxel-scale data).
  if (k < min(dim(W)) && requireNamespace("RSpectra", quietly = TRUE)) {
    sv <- RSpectra::svds(W, k = k, nu = 0, nv = 0)
    d  <- sort(as.numeric(sv$d), decreasing = TRUE)
    return(list(d = d, u = NULL, v = NULL))
  }
  sv <- svd(as.matrix(W), nu = 0, nv = 0)
  d  <- if (k < length(sv$d)) sv$d[seq_len(k)] else sv$d
  list(d = d, u = NULL, v = NULL)
}

#' Principal angles between two subspaces
#'
#' @description
#' Computes principal angles between the column spaces spanned by `U` and `V`.
#'
#' @param U,V Matrices with the same number of rows.
#'
#' @return Numeric vector of angles in radians.
#'
#' @family subspace analysis
#' @seealso [procrustes_distance()] for an alternative subspace distance metric
#' 
#' @examples
#' \donttest{
#' set.seed(1)
#' U <- matrix(rnorm(50), nrow = 10)
#' V <- matrix(rnorm(50), nrow = 10)
#' theta <- subspace_principal_angles(U, V)
#' }
#'
#' @export
subspace_principal_angles <- function(U, V) {
  stopifnot(nrow(U) == nrow(V))
  Uq <- qr.Q(qr(U)); Vq <- qr.Q(qr(V))
  s <- svd(Matrix::t(Uq) %*% Vq)$d
  s[s > 1] <- 1; s[s < 0] <- 0
  acos(s)
}

#' Procrustes distance between two subspaces
#'
#' @description
#' Computes the orthogonal Procrustes distance between subspaces spanned by `U`
#' and `V` using singular values of `U'V`.
#'
#' @param U,V Basis matrices for the two subspaces.
#'
#' @return Nonnegative scalar distance.
#'
#' @family subspace analysis
#' @seealso [subspace_principal_angles()] for computing angles between subspaces
#' 
#' @examples
#' \donttest{
#' set.seed(1)
#' U <- matrix(rnorm(50), nrow = 10)
#' V <- matrix(rnorm(50), nrow = 10)
#' d <- procrustes_distance(U, V)
#' }
#'
#' @export
procrustes_distance <- function(U, V) {
  Uq <- qr.Q(qr(U)); Vq <- qr.Q(qr(V))
  sv <- svd(Matrix::t(Uq) %*% Vq)$d
  k <- ncol(Uq)
  sqrt(max(0, 2*k - 2*sum(sv)))
}

#' Gavish–Donoho-style rank heuristic on whitened singular values
#'
#' @description
#' Estimates noise level and selects `k` by thresholding the singular values `s`
#' of the whitened matrix at `tau = sigma_hat (sqrt(m) + sqrt(n))`.
#'
#' @param W Whitened matrix (T x V) or any numeric matrix.
#'
#' @return A list with `k`, `threshold`, `sigma_hat`, and the spectrum `s`.
#'
#' @family rank selection
#' @seealso [choose_rank_pa()] for parallel analysis approach,
#'   [blocked_cv_recon_error()] for cross-validation approach
#' 
#' @examples
#' \donttest{
#' set.seed(1)
#' W <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' gd <- choose_rank_gd(W)
#' }
#'
#' @export
choose_rank_gd <- function(W) {
  m <- nrow(W); n <- ncol(W)
  colsd <- sqrt(matrixStats::colVars(as.matrix(W)))
  sigma_hat <- stats::median(colsd, na.rm = TRUE) / sqrt(max(m, 1))
  tau <- sigma_hat * (sqrt(m) + sqrt(n))
  s <- svd(W, nu = 0, nv = 0)$d
  k <- sum(s > tau)
  list(k = as.integer(k), threshold = tau, sigma_hat = sigma_hat, s = s)
}

#' Parallel Analysis via time-wise permutation
#'
#' @description
#' Computes a null distribution for the leading singular value by permuting time
#' indices independently within each column and selecting `k` where observed
#' singular values exceed the `(1 - alpha)` quantile of the null.
#'
#' @param W Whitened matrix (T x V).
#' @param B Integer number of permutations (default 100).
#' @param alpha Significance level (default 0.05).
#'
#' @return A list with `k`, `threshold` (null quantile), `null_s1`, and `s`.
#'
#' @family rank selection
#' @seealso [choose_rank_gd()] for Gavish-Donoho approach,
#'   [blocked_cv_recon_error()] for cross-validation approach
#' 
#' @examples
#' \donttest{
#' set.seed(1)
#' W <- matrix(rnorm(200), nrow = 20, ncol = 10)
#' pa <- choose_rank_pa(W, B = 20, alpha = 0.05)
#' }
#'
#' @export
choose_rank_pa <- function(W, B = 100L, alpha = 0.05) {
  W <- as.matrix(W)
  m <- nrow(W); n <- ncol(W)
  s  <- svd(W, nu = 0, nv = 0)$d
  s1_null <- numeric(B)
  for (b in seq_len(B)) {
    Wb <- W
    for (j in seq_len(n)) {
      Wb[, j] <- Wb[sample.int(m, size = m, replace = FALSE), j]
    }
    s1_null[b] <- svd(Wb, nu = 0, nv = 0)$d[1]
  }
  thr <- stats::quantile(s1_null, probs = 1 - alpha, type = 7, names = FALSE)
  k <- sum(s > thr)
  list(k = as.integer(k), threshold = thr, null_s1 = s1_null, s = s)
}

#' Blocked cross-validated reconstruction error under genpca geometry
#'
#' @description
#' Evaluates reconstruction error by holding out contiguous time blocks and
#' projecting hold-out data onto the `k`-dimensional column space learned on the
#' training block, all under the whitened geometry.
#'
#' @param X,A,M Data and metrics as in `whitened_matrix()`.
#' @param k Number of components.
#' @param nfold Number of folds (hold-out blocks).
#' @param block_frac Fraction of time points per hold-out block.
#'
#' @return Data frame with per-fold errors and the mean.
#'
#' @family rank selection
#' @seealso [choose_rank_gd()] for Gavish-Donoho approach,
#'   [choose_rank_pa()] for parallel analysis approach
#' 
#' @examples
#' \donttest{
#' set.seed(1)
#' X <- matrix(rnorm(500), nrow = 50, ncol = 10)
#' A <- Matrix::Diagonal(ncol(X))
#' M <- Matrix::Diagonal(nrow(X))
#' cv <- blocked_cv_recon_error(X, A, M, k = 3)
#' }
#'
#' @export
blocked_cv_recon_error <- function(X, A, M, k, nfold = 5L, block_frac = 0.2) {
  Tlen <- nrow(X)
  # Disjoint, contiguous hold-out blocks: stride >= block size guarantees folds
  # never overlap.
  stride <- max(1L, floor(Tlen / nfold))
  bsz    <- max(1L, min(floor(Tlen * block_frac), stride))
  starts <- 1L + (seq_len(nfold) - 1L) * stride
  errs   <- numeric(nfold)

  RA <- chol(A)

  for (f in seq_len(nfold)) {
    s0 <- starts[f]
    idx_hold  <- seq.int(s0, length.out = min(bsz, Tlen - s0 + 1L))
    idx_train <- setdiff(seq_len(Tlen), idx_hold)

    # Whiten rows with the Cholesky of the relevant *submatrix* of M; a
    # submatrix of chol(M) is not the Cholesky of the submatrix of M.
    RM_tr <- chol(Matrix::forceSymmetric(M[idx_train, idx_train, drop = FALSE]))
    RM_ho <- chol(Matrix::forceSymmetric(M[idx_hold,  idx_hold,  drop = FALSE]))

    W_tr <- (Matrix::t(RM_tr)) %*% X[idx_train, , drop = FALSE] %*% RA
    sv   <- svd(as.matrix(W_tr), nu = k, nv = k)
    Vk   <- sv$v[, seq_len(k), drop = FALSE]

    W_ho  <- (Matrix::t(RM_ho)) %*% X[idx_hold, , drop = FALSE] %*% RA
    # Project hold-out onto the column subspace learned from training.
    W_rec <- W_ho %*% Vk %*% t(Vk)

    errs[f] <- sqrt(sum(as.matrix(W_ho - W_rec)^2) / length(W_ho))
  }
  data.frame(fold = seq_len(nfold), error = errs, mean_error = mean(errs))
}
