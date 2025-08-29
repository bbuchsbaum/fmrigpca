pkgname <- "fmrigpca"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('fmrigpca')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("aggregate_tissue_to_parcels")
### * aggregate_tissue_to_parcels

flush(stderr()); flush(stdout())

### Name: aggregate_tissue_to_parcels
### Title: Aggregate GM/WM/CSF tissue probabilities to parcels
### Aliases: aggregate_tissue_to_parcels

### ** Examples

## Not run: 
##D agg <- aggregate_tissue_to_parcels(gm, wm, csf, parc)
## End(Not run)




cleanEx()
nameEx("blocked_cv_recon_error")
### * blocked_cv_recon_error

flush(stderr()); flush(stdout())

### Name: blocked_cv_recon_error
### Title: Blocked cross-validated reconstruction error under genpca
###   geometry
### Aliases: blocked_cv_recon_error

### ** Examples

## Not run: 
##D cv <- blocked_cv_recon_error(X, A, M, k = 5)
## End(Not run)




cleanEx()
nameEx("build_temporal_metric")
### * build_temporal_metric

flush(stderr()); flush(stdout())

### Name: build_temporal_metric
### Title: Build temporal covariance metric for generalized PCA
### Aliases: build_temporal_metric

### ** Examples

## Not run: 
##D # Basic usage with AR(1) whitening
##D M <- build_temporal_metric(nv_run, p = 1L)
##D 
##D # With motion correction
##D M <- build_temporal_metric(nv_run, wm_mask = wm_vol, 
##D                            FD = fd_values, DVARS = dvars_values,
##D                            lambda_t = 0.5)
##D 
##D # No AR modeling for task fMRI with explicit task regressors
##D M <- build_temporal_metric(nv_run, p = 0, lambda_t = 0.3)
## End(Not run)




cleanEx()
nameEx("build_temporal_metric_parcel")
### * build_temporal_metric_parcel

flush(stderr()); flush(stdout())

### Name: build_temporal_metric_parcel
### Title: Build temporal covariance metric for parcellated fMRI data
### Aliases: build_temporal_metric_parcel

### ** Examples

## Not run: 
##D # Build temporal metric with AR(1) whitening and motion regressors
##D M <- build_temporal_metric_parcel(cnv_run, wm_parcels = c(1,2,3),
##D                                   FD = fd_vals, DVARS = dvars_vals)
## End(Not run)




cleanEx()
nameEx("choose_rank_gd")
### * choose_rank_gd

flush(stderr()); flush(stdout())

### Name: choose_rank_gd
### Title: Gavish–Donoho-style rank heuristic on whitened singular values
### Aliases: choose_rank_gd

### ** Examples

## Not run: 
##D gd <- choose_rank_gd(W)
## End(Not run)




cleanEx()
nameEx("choose_rank_pa")
### * choose_rank_pa

flush(stderr()); flush(stdout())

### Name: choose_rank_pa
### Title: Parallel Analysis via time-wise permutation
### Aliases: choose_rank_pa

### ** Examples

## Not run: 
##D pa <- choose_rank_pa(W, B = 200, alpha = 0.05)
## End(Not run)




cleanEx()
nameEx("compute_tsnr")
### * compute_tsnr

flush(stderr()); flush(stdout())

### Name: compute_tsnr
### Title: Temporal SNR (tSNR) per voxel across runs
### Aliases: compute_tsnr

### ** Examples

## Not run: 
##D tsnr <- compute_tsnr(list(nv1, nv2))
## End(Not run)




cleanEx()
nameEx("compute_tsnr_parcel")
### * compute_tsnr_parcel

flush(stderr()); flush(stdout())

### Name: compute_tsnr_parcel
### Title: Parcel-level tSNR across runs (ClusteredNeuroVec)
### Aliases: compute_tsnr_parcel

### ** Examples

## Not run: 
##D # Compute parcel tSNR across multiple runs
##D tsnr_p <- compute_tsnr_parcel(cnv_list)
## End(Not run)




cleanEx()
nameEx("estimate_ar_whitener")
### * estimate_ar_whitener

flush(stderr()); flush(stdout())

### Name: estimate_ar_whitener
### Title: Estimate a design-free AR(p) whitener from a run
### Aliases: estimate_ar_whitener

### ** Examples

## Not run: 
##D nv_run <- neuroim2::simulate_fmri(mask, n_time = 50, seed = 1)
##D arw   <- estimate_ar_whitener(nv_run, p = 1L)
## End(Not run)




cleanEx()
nameEx("estimate_ar_whitener_parcel")
### * estimate_ar_whitener_parcel

flush(stderr()); flush(stdout())

### Name: estimate_ar_whitener_parcel
### Title: Estimate AR(p) whitener from parcel-level data
### Aliases: estimate_ar_whitener_parcel

### ** Examples

## Not run: 
##D # Estimate AR(1) whitener using WM parcels
##D arw <- estimate_ar_whitener_parcel(cnv_run, wm_parcels = c(1,2,3), p = 1)
## End(Not run)




cleanEx()
nameEx("fit_subject_genpca")
### * fit_subject_genpca

flush(stderr()); flush(stdout())

### Name: fit_subject_genpca
### Title: Fit generalized PCA across multiple runs (voxel level)
### Aliases: fit_subject_genpca

### ** Examples

## Not run: 
##D fit <- fit_subject_genpca(nv_list, mask, gm, wm, csf, k = 5)
## End(Not run)




cleanEx()
nameEx("fit_subject_genpca_parcel")
### * fit_subject_genpca_parcel

flush(stderr()); flush(stdout())

### Name: fit_subject_genpca_parcel
### Title: Fit generalized PCA at the parcel level (ClusteredNeuroVec)
### Aliases: fit_subject_genpca_parcel

### ** Examples

## Not run: 
##D # Fit genpca on parcellated data with tissue weights
##D fit <- fit_subject_genpca_parcel(cnv_list, parc_vol = parcellation,
##D                                  gm_vol = gm, wm_vol = wm, csf_vol = csf,
##D                                  k = 20)
## End(Not run)




cleanEx()
nameEx("fit_subject_metapca")
### * fit_subject_metapca

flush(stderr()); flush(stdout())

### Name: fit_subject_metapca
### Title: Combine per-run genpca fits via meta-PCA (MFA)
### Aliases: fit_subject_metapca

### ** Examples

## Not run: 
##D # Per-run fits (each list contains a single run)
##D fits <- lapply(seq_along(nv_list), function(r) {
##D   fit_subject_genpca(
##D     nv_list = list(nv_list[[r]]),
##D     mask_vol = mask,
##D     gm_vol = gm, wm_vol = wm, csf_vol = csf,
##D     k = 5, p_ar = 1L
##D   )
##D })
##D # Combine via meta-PCA (MFA weighting)
##D meta <- fit_subject_metapca(fits, k = 3)
##D str(meta$v)  # consensus loadings (voxels x k)
## End(Not run)




cleanEx()
nameEx("make_A_from_parcels")
### * make_A_from_parcels

flush(stderr()); flush(stdout())

### Name: make_A_from_parcels
### Title: Build parcel-level column metric A from aggregated tissues and
###   Laplacian
### Aliases: make_A_from_parcels

### ** Examples

## Not run: 
##D # Build parcel column metric with tissue weights and spatial smoothing
##D A <- make_A_from_parcels(gm_p, wm_p, csf_p, Lp, tsnr_p = tsnr_vals,
##D                         lambda_s = 0.5)
## End(Not run)




cleanEx()
nameEx("make_A_from_vols")
### * make_A_from_vols

flush(stderr()); flush(stdout())

### Name: make_A_from_vols
### Title: Build voxel (column) metric A from tissue maps, Laplacian, and
###   tSNR
### Aliases: make_A_from_vols

### ** Examples

## Not run: 
##D A <- make_A_from_vols(gm, wm, csf, L, tsnr = ts, mask_idx = idx)
## End(Not run)




cleanEx()
nameEx("make_frame_weights")
### * make_frame_weights

flush(stderr()); flush(stdout())

### Name: make_frame_weights
### Title: Frame-wise robust weights from FD and DVARS
### Aliases: make_frame_weights

### ** Examples

## Not run: 
##D W <- make_frame_weights(FD = runif(30, 0, 0.5), DVARS = runif(30, 0.5, 1.5))
## End(Not run)




cleanEx()
nameEx("make_laplacian")
### * make_laplacian

flush(stderr()); flush(stdout())

### Name: make_laplacian
### Title: Build a spatial graph Laplacian from a NeuroVol mask
### Aliases: make_laplacian

### ** Examples

## Not run: 
##D dims <- c(6,6,4)
##D arr  <- array(FALSE, dims); arr[3:4,3:4,2:3] <- TRUE
##D mv   <- neuroim2::NeuroVol(arr, neuroim2::NeuroSpace(dims, c(2,2,2)))
##D L    <- make_laplacian(mv, k = 4, sigma = 2, normalized = TRUE)
##D Matrix::isSymmetric(L)
## End(Not run)




cleanEx()
nameEx("make_parcel_laplacian")
### * make_parcel_laplacian

flush(stderr()); flush(stdout())

### Name: make_parcel_laplacian
### Title: Build a parcel-level Laplacian
### Aliases: make_parcel_laplacian

### ** Examples

## Not run: 
##D Lp <- make_parcel_laplacian(parc_vol = parc, k = 6, sigma = 2)
## End(Not run)




cleanEx()
nameEx("make_temporal_penalty")
### * make_temporal_penalty

flush(stderr()); flush(stdout())

### Name: make_temporal_penalty
### Title: Second-difference temporal smoothing penalty
### Aliases: make_temporal_penalty

### ** Examples

## Not run: 
##D H <- make_temporal_penalty(30, lambda_t = 0.5)
## End(Not run)




cleanEx()
nameEx("parcel_centroids_from_labels")
### * parcel_centroids_from_labels

flush(stderr()); flush(stdout())

### Name: parcel_centroids_from_labels
### Title: Parcel centroids from a label volume
### Aliases: parcel_centroids_from_labels

### ** Examples

## Not run: 
##D pc <- parcel_centroids_from_labels(parc_vol)
## End(Not run)




cleanEx()
nameEx("procrustes_distance")
### * procrustes_distance

flush(stderr()); flush(stdout())

### Name: procrustes_distance
### Title: Procrustes distance between two subspaces
### Aliases: procrustes_distance

### ** Examples

## Not run: 
##D d <- procrustes_distance(U, V)
## End(Not run)




cleanEx()
nameEx("subspace_principal_angles")
### * subspace_principal_angles

flush(stderr()); flush(stdout())

### Name: subspace_principal_angles
### Title: Principal angles between two subspaces
### Aliases: subspace_principal_angles

### ** Examples

## Not run: 
##D theta <- subspace_principal_angles(U, V)
## End(Not run)




cleanEx()
nameEx("whitened_matrix")
### * whitened_matrix

flush(stderr()); flush(stdout())

### Name: whitened_matrix
### Title: Form a whitened matrix for genpca geometry
### Aliases: whitened_matrix

### ** Examples

## Not run: 
##D W <- whitened_matrix(X, A, M)
## End(Not run)




cleanEx()
nameEx("whitened_svd")
### * whitened_svd

flush(stderr()); flush(stdout())

### Name: whitened_svd
### Title: Truncated SVD of the whitened matrix
### Aliases: whitened_svd

### ** Examples

## Not run: 
##D ws <- whitened_svd(X, A, M, k = 10)
## End(Not run)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
