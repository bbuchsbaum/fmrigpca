#' Temporal SNR (tSNR) per voxel across runs
#'
#' @description
#' Computes per-voxel temporal signal-to-noise ratio across a list of `NeuroVec`
#' runs by aggregating first and second moments over time and runs.
#'
#' @details
#' For each voxel, the mean and variance are computed over concatenated time
#' points across runs. tSNR is `max(1e-6, mean) / max(1e-6, sd)`, then lower-
#' bounded at `1e-6` to avoid degenerate weights.
#'
#' @param nv_list List of `neuroim2::NeuroVec` runs.
#'
#' @return Numeric vector of length V (voxels) with tSNR values.
#'
#' @family spatial metrics
#' @seealso [compute_tsnr_parcel()] for parcel-level tSNR computation
#'
#' @examples
#' \donttest{
#' # Create simple test data
#' dims <- c(10, 10, 5, 20)  # 4D: x, y, z, time
#' space <- neuroim2::NeuroSpace(dims,
#'                                spacing = c(1, 1, 1),
#'                                origin = c(0, 0, 0))
#' mask <- array(TRUE, dims[1:3])  # 3D mask
#'
#' # Create two runs with random data
#' data1 <- array(rnorm(prod(dims)), dims)
#' data2 <- array(rnorm(prod(dims)), dims)
#' nv1 <- neuroim2::NeuroVec(data1, space, mask = mask)
#' nv2 <- neuroim2::NeuroVec(data2, space, mask = mask)
#' 
#' tsnr <- compute_tsnr(list(nv1, nv2))
#' }
#'
#' @export
compute_tsnr <- function(nv_list) {
  stopifnot(length(nv_list) >= 1L)
  X0 <- as.matrix(nv_list[[1L]])
  V  <- nrow(X0)
  s1 <- numeric(V)
  s2 <- numeric(V)
  Ttot <- 0L
  for (nv in nv_list) {
    X <- as.matrix(nv)   # V x T
    s1 <- s1 + Matrix::rowSums(X)
    s2 <- s2 + Matrix::rowSums(X^2)
    Ttot <- Ttot + ncol(X)
  }
  mu <- s1 / Ttot
  var <- pmax(0, s2 / Ttot - mu^2)
  sdv <- sqrt(var)
  pmax(1e-6, pmax(1e-6, mu) / pmax(1e-6, sdv))
}

#' Build spatial column metric A from tissue maps, Laplacian, and tSNR
#'
#' @description
#' Constructs a positive definite spatial column metric `A` that encodes spatial
#' smoothness (graph Laplacian) and voxelwise weights derived from tissue
#' probabilities and temporal SNR. Adds a small ridge for numerical stability.
#'
#' @details
#' Let `w = gm^alpha * tsnr^beta * (wm + csf)^{-gamma}` on in-mask voxels. Define
#' `S = I + lambda_s * L` where `L` is a (normalized) Laplacian. Then
#' `A = diag(sqrt(w)) S diag(sqrt(w)) + tau * I`. The `mask_idx` selects in-mask
#' voxel indices and must align with `nrow(L)`. If `tsnr` is `NULL`, it defaults
#' to 1. Tissue maps must share space with the mask.
#'
#' @param gm_vol,wm_vol,csf_vol `NeuroVol` tissue probability maps.
#' @param L Laplacian `Matrix` over in-mask voxels.
#' @param tsnr Optional numeric tSNR vector over in-mask voxels.
#' @param mask_idx Integer indices of in-mask voxels (length equals `nrow(L)`).
#' @param alpha,beta,gamma Exponents for GM, tSNR, and (WM+CSF) weighting.
#' @param lambda_s Spatial regularization weight multiplying `L`.
#' @param tau Small ridge added to the diagonal (default 1e-6).
#'
#' @return Symmetric positive definite `Matrix` of size `V_in x V_in`.
#'
#' @family spatial metrics
#' @seealso [build_spatial_metric_parcel()] for parcel-level column metric construction
#'
#' @examples
#' \donttest{
#' # Create small example volumes
#' dims <- c(5, 5, 3)
#' space <- neuroim2::NeuroSpace(dims, c(1, 1, 1))
#' 
#' # Create tissue probability maps
#' gm_data <- array(runif(prod(dims), 0.5, 1), dims)
#' wm_data <- array(runif(prod(dims), 0, 0.3), dims)
#' csf_data <- array(runif(prod(dims), 0, 0.2), dims)
#' 
#' gm <- neuroim2::NeuroVol(gm_data, space)
#' wm <- neuroim2::NeuroVol(wm_data, space)
#' csf <- neuroim2::NeuroVol(csf_data, space)
#' 
#' # Create mask and Laplacian
#' mask <- neuroim2::NeuroVol(gm_data > 0.3, space)
#' L <- make_laplacian(mask, k = 6)
#' mask_idx <- which(neuroim2::values(mask) > 0)
#' 
#' A <- build_spatial_metric(gm, wm, csf, L, mask_idx = mask_idx)
#' }
#'
#' @export
build_spatial_metric <- function(gm_vol, wm_vol, csf_vol, L,
                             tsnr = NULL, mask_idx = NULL,
                             alpha = 1.0, beta = 0.5, gamma = 1.0,
                             lambda_s = 0.5, tau = 1e-6) {
  stopifnot(inherits(gm_vol, "NeuroVol"),
            inherits(wm_vol, "NeuroVol"),
            inherits(csf_vol, "NeuroVol"))
  sp <- neuroim2::space(gm_vol)
  if (!identical(sp, neuroim2::space(wm_vol)) || !identical(sp, neuroim2::space(csf_vol)))
    stop("Tissue maps must share the same NeuroSpace.")

  gm_all  <- as.numeric(neuroim2::values(gm_vol))
  wm_all  <- as.numeric(neuroim2::values(wm_vol))
  csf_all <- as.numeric(neuroim2::values(csf_vol))

  if (is.null(mask_idx)) {
    mask_idx <- which(gm_all + wm_all + csf_all > 0)
    if (length(mask_idx) != nrow(L)) {
      stop("mask_idx must align with the Laplacian L dimension.")
    }
  }

  gm  <- gm_all[mask_idx]
  wm  <- wm_all[mask_idx]
  csf <- csf_all[mask_idx]
  if (is.null(tsnr)) {
    tsnr <- rep(1, length(gm))
  } else if (length(tsnr) != length(gm)) {
    tsnr <- tsnr[mask_idx]
  }

  wv <- (pmax(gm, 1e-6)^alpha) * (pmax(tsnr, 1e-6)^beta) * ((pmax(wm + csf, 1e-6))^(-gamma))
  S  <- Matrix::Diagonal(n = length(wv)) + lambda_s * L
  A  <- Matrix::forceSymmetric((Matrix::Diagonal(x = sqrt(wv)) %*% S %*% Matrix::Diagonal(x = sqrt(wv))))
  A + Matrix::Diagonal(nrow(A)) * tau
}
