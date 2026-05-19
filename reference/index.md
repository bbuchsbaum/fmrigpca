# Package index

## All functions

- [`aggregate_tissue_to_parcels()`](https://bbuchsbaum.github.io/fmrigpca/reference/aggregate_tissue_to_parcels.md)
  : Aggregate GM/WM/CSF tissue probabilities to parcels
- [`blocked_cv_recon_error()`](https://bbuchsbaum.github.io/fmrigpca/reference/blocked_cv_recon_error.md)
  : Blocked cross-validated reconstruction error under genpca geometry
- [`build_spatial_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric.md)
  : Build spatial column metric A from tissue maps, Laplacian, and tSNR
- [`build_spatial_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric_parcel.md)
  : Build parcel-level spatial column metric A from aggregated tissues
  and Laplacian
- [`build_temporal_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric.md)
  : Build temporal covariance metric for generalized PCA
- [`build_temporal_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric_parcel.md)
  : Build temporal covariance metric for parcellated fMRI data
- [`choose_rank_gd()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_gd.md)
  : Gavish–Donoho-style rank heuristic on whitened singular values
- [`choose_rank_pa()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_pa.md)
  : Parallel Analysis via time-wise permutation
- [`compute_tsnr()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr.md)
  : Temporal SNR (tSNR) per voxel across runs
- [`compute_tsnr_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr_parcel.md)
  : Parcel-level tSNR across runs (ClusteredNeuroVec)
- [`estimate_ar_whitener()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener.md)
  : Estimate a design-free AR(p) whitener from a run
- [`estimate_ar_whitener_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener_parcel.md)
  : Estimate AR(p) whitener from parcel-level data
- [`fit_subject_genpca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca.md)
  : Fit generalized PCA across multiple runs (voxel level)
- [`fit_subject_genpca_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca_parcel.md)
  : Fit generalized PCA at the parcel level (ClusteredNeuroVec)
- [`fit_subject_metapca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_metapca.md)
  : Combine per-run genpca fits via meta-PCA (MFA)
- [`make_frame_weights()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md)
  : Frame-wise robust weights from FD and DVARS
- [`make_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_laplacian.md)
  : Build a spatial graph Laplacian from a NeuroVol mask
- [`make_parcel_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_parcel_laplacian.md)
  : Build a parcel-level Laplacian
- [`make_temporal_penalty()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_temporal_penalty.md)
  : Second-difference temporal smoothing penalty
- [`parcel_centroids_from_labels()`](https://bbuchsbaum.github.io/fmrigpca/reference/parcel_centroids_from_labels.md)
  : Parcel centroids from a label volume
- [`procrustes_distance()`](https://bbuchsbaum.github.io/fmrigpca/reference/procrustes_distance.md)
  : Procrustes distance between two subspaces
- [`subspace_principal_angles()`](https://bbuchsbaum.github.io/fmrigpca/reference/subspace_principal_angles.md)
  : Principal angles between two subspaces
- [`whitened_matrix()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_matrix.md)
  : Form a whitened matrix for genpca geometry
- [`whitened_svd()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_svd.md)
  : Truncated SVD of the whitened matrix
