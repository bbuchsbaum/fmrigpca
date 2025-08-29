#' Build a spatial graph Laplacian from a NeuroVol mask
#'
#' @description
#' Constructs a sparse, symmetric voxel–wise graph Laplacian for a masked volume
#' using k–nearest neighbors with a Gaussian kernel. Optionally returns the
#' normalized Laplacian and scales by the largest eigenvalue for numerical
#' stability.
#'
#' @details
#' Voxels inside the mask are converted to integer coordinates (i,j,k), then a
#' k–NN graph is built with edge weights `exp(-d^2 / (2*sigma^2))` using the
#' `neighborweights` package (functions `knn_weights` or `spatial_adjacency`).
#' The unnormalized Laplacian is `L = D - W`.
#' The normalized variant is `L = I - D^{-1/2} W D^{-1/2}`. In both cases, `L`
#' is divided by its largest eigenvalue when positive and finite, keeping the
#' spectrum within a stable range (e.g., 0 to 2 for the normalized Laplacian).
#'
#' Edge cases are handled: a single-voxel mask yields a 1x1 identity; empty
#' masks are rejected. If `k` is too large, it is clipped to `n-1` neighbors.
#'
#' @param mask_vol A `neuroim2::NeuroVol` logical/numeric mask defining in-mask voxels.
#' @param k Integer; number of neighbors in the k–NN graph (default 6). Common
#'   choices: 6 for face connectivity, 18 for face+edge, 26 for full 3D.
#' @param sigma Numeric; Gaussian kernel bandwidth for edge weights (default 2.5).
#'   Smaller values create sparser graphs; larger values smooth more broadly.
#' @param normalized Logical; return the normalized Laplacian if `TRUE` (default).
#'   Normalized version has better numerical properties and bounded spectrum.
#'
#' @return A sparse, symmetric `Matrix` of size `n_voxels x n_voxels` representing
#'   the (optionally normalized) graph Laplacian over in-mask voxels.
#'
#' @family spatial metrics
#' @seealso [make_parcel_laplacian()] for parcel-level Laplacian construction
#' 
#' @examples
#' \dontrun{
#' dims <- c(6,6,4)
#' arr  <- array(FALSE, dims); arr[3:4,3:4,2:3] <- TRUE
#' mv   <- neuroim2::NeuroVol(arr, neuroim2::NeuroSpace(dims, c(2,2,2)))
#' L    <- make_laplacian(mv, k = 4, sigma = 2, normalized = TRUE)
#' Matrix::isSymmetric(L)
#' }
#'
#' @export
make_laplacian <- function(mask_vol, k = 6, sigma = 2.5,
                           normalized = TRUE) {
  stopifnot(inherits(mask_vol, "NeuroVol"))
  vals <- neuroim2::values(mask_vol)
  idx  <- which(vals != 0)
  dims <- dim(mask_vol)
  if (length(idx) == 0L) stop("Mask has no nonzero voxels.")

  Vx <- dims[1]; Vy <- dims[2]; Vz <- if (length(dims) >= 3) dims[3] else 1
  lin <- idx - 1L
  ii <- (lin %% Vx) + 1L
  jj <- ((lin %/% Vx) %% Vy) + 1L
  kk <- (lin %/% (Vx * Vy)) + 1L
  coords <- cbind(ii, jj, kk)

  # Use neighborweights for efficient k-NN graph construction if available
  W <- if (requireNamespace("neighborweights", quietly = TRUE)) {
    tryCatch({
      g <- neighborweights::graph_weights(
        X = coords,
        k = k,
        neighbor_mode = "knn",
        weight_mode = "heat",
        type = "normal",
        sigma = sigma
      )
      # Extract adjacency matrix from neighbor_graph object
      neighborweights::adjacency(g)
    }, error = function(e) {
      # Fallback to internal implementation if neighborweights fails
      .knn_gaussian(coords, k = k, sigma = sigma, symmetric = TRUE)
    })
  } else {
    .knn_gaussian(coords, k = k, sigma = sigma, symmetric = TRUE)
  }
  W <- Matrix::forceSymmetric(W)

  # Handle special case of single node
  if (nrow(W) == 1) {
    L <- Matrix::Matrix(1, nrow = 1, ncol = 1)
    return(L)
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

  lam_max <- if (requireNamespace("RSpectra", quietly = TRUE)) {
    # Convert to regular matrix if very small, otherwise use RSpectra
    if (nrow(L) < 100) {
      max(Re(eigen(as.matrix(L), only.values = TRUE)$values))
    } else {
      as.numeric(RSpectra::eigs_sym(L, 1, which = "LM")$values)
    }
  } else {
    max(Re(eigen(as.matrix(L), only.values = TRUE)$values))
  }
  if (is.finite(lam_max) && lam_max > 0) L <- L / lam_max
  L
}

.knn_gaussian <- function(coords, k = 6, sigma = 2.5, symmetric = TRUE) {
  n <- nrow(coords)
  
  # Handle edge case where k is too large for the number of points
  if (n == 1) {  # Single point case
    W <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(1, 1))
    return(W)
  }
  
  if (k >= n) {
    k <- n - 1
  }
  
  if (k == 0) {  # Only happens if n == 1, but already handled above
    W <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
    return(W)
  }
  
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn  <- RANN::nn2(coords, query = coords, k = min(k + 1L, n))
    idx <- nn$nn.idx[, -1, drop = FALSE]
    dst <- nn$nn.dists[, -1, drop = FALSE]
  } else {
    D <- as.matrix(dist(coords))
    diag(D) <- Inf
    idx <- t(apply(D, 1L, function(row) order(row)[1:k]))
    dst <- matrix(NA_real_, nrow = n, ncol = k)
    for (i in seq_len(n)) dst[i, ] <- D[i, idx[i, ]]
  }
  i <- rep.int(seq_len(n), times = k)
  j <- as.vector(idx)
  w <- exp(-(as.vector(dst)^2) / (2 * sigma^2))
  W <- Matrix::sparseMatrix(i = i, j = j, x = w, dims = c(n, n))
  if (symmetric) W <- pmax(W, Matrix::t(W))
  W
}
