#' Fit generalized PCA across multiple runs (voxel level)
#'
#' @description
#' Builds a voxel-level column metric `A` and per-run row metrics `M_r`, stacks
#' runs in time, and fits `genpca` with `A` and block-diagonal `M`. Adds
#' `scores` (projected time scores) for convenience.
#'
#' @details
#' Validates that all `NeuroVec` runs share the same spatial `NeuroSpace` as the
#' mask and tissue maps. Constructs a Laplacian `L` via `make_laplacian`, a
#' tissue/tSNR-weighted `A` via `build_spatial_metric`, and per-run `M_r` via
#' `make_M_for_run` using optional FD/DVARS. Rows (time) are stacked across runs
#' and `genpca::genpca` is called with `ncomp = k`. The result gains a `scores`
#' component via `multivarious::scores`.
#'
#' @param nv_list List of `NeuroVec` runs.
#' @param mask_vol `NeuroVol` binary mask.
#' @param gm_vol,wm_vol,csf_vol `NeuroVol` tissue probability maps (same space as mask).
#' @param FD_list,DVARS_list Optional lists of per-run motion metrics.
#'   Each list should contain numeric vectors with one element per list 
#'   corresponding to each run in nv_list. FD is framewise displacement (mm)
#'   and DVARS is the temporal derivative of signal variance. These metrics
#'   are used to down-weight high-motion frames. See \code{\link{make_frame_weights}}
#'   for detailed definitions.
#' @param k Number of components to fit (default 20). Start with modest values
#'   (5-10) for exploration, increase based on rank selection diagnostics.
#' @param p_ar AR order for row whitener (default 1). Use 0 to skip AR modeling,
#'   higher values (2-3) for data with stronger temporal autocorrelation.
#' @param lambda_s Spatial regularization weight (default 0.5). Controls smoothness
#'   via the graph Laplacian; higher values enforce more spatial coherence.
#' @param lambda_t Temporal regularization weight (default 0.3). Controls temporal
#'   smoothness; higher values penalize rapid fluctuations.
#' @param lap_k Number of nearest neighbors for Laplacian construction (default 6).
#'   Typically 6-8 for 3D connectivity; higher values increase spatial coupling.
#' @param lap_sigma Gaussian kernel bandwidth for edge weights (default 2.5).
#'   Controls the decay of spatial influence with distance.
#'
#' @return A `genpca` object with loadings `v` and added `scores` (stacked time x k).
#'
#' @family main fitting
#' @seealso [fit_subject_genpca_parcel()] for parcel-level analysis,
#'   [fit_subject_metapca()] for combining multiple fits
#' 
#' @examples
#' \dontrun{
#' fit <- fit_subject_genpca(nv_list, mask, gm, wm, csf, k = 5)
#' }
#'
#' @export
fit_subject_genpca <- function(nv_list, mask_vol, gm_vol, wm_vol, csf_vol,
                               FD_list = NULL, DVARS_list = NULL,
                               k = 20, p_ar = 1L, lambda_s = 0.5, lambda_t = 0.3,
                               lap_k = 6, lap_sigma = 2.5) {
  stopifnot(length(nv_list) >= 1L, inherits(mask_vol, "NeuroVol"))
  sp <- neuroim2::space(mask_vol)
  # Check that the first 3 dimensions match (NeuroVec is 4D, mask is 3D)
  sp_check <- function(nv) {
    nv_sp <- neuroim2::space(nv)
    # Compare first 3 dimensions
    identical(nv_sp@dim[1:3], sp@dim) && 
    identical(nv_sp@spacing[1:3], sp@spacing) &&
    identical(nv_sp@origin[1:3], sp@origin)
  }
  if (!all(vapply(nv_list, sp_check, logical(1L))))
    stop("All runs must share the same spatial dimensions as mask_vol.")
  if (!identical(sp, neuroim2::space(gm_vol)) ||
      !identical(sp, neuroim2::space(wm_vol)) ||
      !identical(sp, neuroim2::space(csf_vol)))
    stop("Tissue maps must share the same NeuroSpace as mask_vol.")

  mask_vals <- neuroim2::values(mask_vol)
  mask_idx  <- which(mask_vals != 0)
  if (length(mask_idx) == 0L) stop("Mask has no voxels.")

  L    <- make_laplacian(mask_vol, k = lap_k, sigma = lap_sigma, normalized = TRUE)
  tsnr_full <- compute_tsnr(nv_list)
  A    <- build_spatial_metric(gm_vol, wm_vol, csf_vol, L, tsnr = tsnr_full, mask_idx = mask_idx,
                           lambda_s = lambda_s)

  X_rows <- list(); M_blks <- list()
  for (r in seq_along(nv_list)) {
    nv <- nv_list[[r]]
    Xr_full <- t(as.matrix(nv))           # T_r x V_full
    Xr <- Xr_full[, mask_idx, drop = FALSE]
    FD_r    <- if (!is.null(FD_list))    FD_list[[r]]    else NULL
    DVARS_r <- if (!is.null(DVARS_list)) DVARS_list[[r]] else NULL
    Mr <- build_temporal_metric(nv, wm_mask = wm_vol, p = p_ar, FD = FD_r, DVARS = DVARS_r,
                                lambda_t = lambda_t, ridge = 1e-6)
    X_rows[[r]] <- Xr
    M_blks[[r]] <- Mr
  }
  X_stack <- do.call(rbind, X_rows)
  M_blk   <- Matrix::bdiag(M_blks)

  result <- genpca::genpca(X_stack, A = A, M = M_blk, ncomp = k, preproc = multivarious::center())
  
  # Add scores component for compatibility
  result$scores <- multivarious::scores(result)
  
  result
}

#' Combine per-run genpca fits via meta-PCA (MFA)
#'
#' @description
#' Aggregates multiple `genpca` fits (e.g., per-run, per-session, or per-subject)
#' into a single consensus component space using `subpca::metapca` with Multiple
#' Factor Analysis (MFA) weighting. This is a “PCA of PCAs” that operates on the
#' loadings from each individual fit, balancing contributions across inputs.
#'
#' @details
#' Two common workflows are supported in fmrigpca:
#' 1) Joint fit across runs: stack runs in time, build a block-diagonal row
#'    metric `M`, and fit once via [fit_subject_genpca()]. This shares a single
#'    spatial column metric `A` and exploits all timepoints jointly.
#' 2) Meta-PCA across individual fits: fit each run (or session/subject) 
#'    separately to obtain `genpca` objects, then combine them with meta-PCA. 
#'
#' Meta-PCA uses MFA weighting to avoid any single fit dominating due to scale
#' or sample size differences, and returns a `metapca` object with combined
#' loadings (accessible as `$v`). Use this when:
#' - Runs/sessions are heterogeneous and you prefer to fit them independently
#'   (e.g., different preprocessing, durations, or acquisition settings)
#' - Combining across subjects after per-subject genpca fits
#'
#' Assumptions and tips:
#' - Loadings should be defined on the same voxel/parcel index set (e.g., same
#'   mask and order). If they differ, align to a common subset before fitting.
#' - Keep `k` small and comparable across inputs; `metapca` selects a 
#'   consensus space up to `ncomp = k`.
#' - See the subpca package for algorithmic details of MFA combination.
#'
#' @param fits List of `genpca` fit objects to combine.
#' @param k Number of meta-components to return (consensus dimensionality).
#'
#' @return A `metapca` object with combined loadings `v` and related metadata.
#'
#' @family main fitting
#' @seealso [fit_subject_genpca()] and [fit_subject_genpca_parcel()] for
#'   generating individual fits to combine
#'
#' @examples
#' \dontrun{
#' # Per-run fits (each list contains a single run)
#' fits <- lapply(seq_along(nv_list), function(r) {
#'   fit_subject_genpca(
#'     nv_list = list(nv_list[[r]]),
#'     mask_vol = mask,
#'     gm_vol = gm, wm_vol = wm, csf_vol = csf,
#'     k = 5, p_ar = 1L
#'   )
#' })
#' # Combine via meta-PCA (MFA weighting)
#' meta <- fit_subject_metapca(fits, k = 3)
#' str(meta$v)  # consensus loadings (voxels x k)
#' }
#'
#' @export
fit_subject_metapca <- function(fits, k = 20) {
  if (!requireNamespace("subpca", quietly = TRUE)) {
    stop("subpca not installed. Install from GitHub: bbuchsbaum/subpca")
  }
  subpca::metapca(fits, ncomp = k, combine = "MFA")
}
