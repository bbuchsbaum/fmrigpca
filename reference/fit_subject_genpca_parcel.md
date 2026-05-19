# Fit generalized PCA at the parcel level (ClusteredNeuroVec)

Builds a parcel Laplacian and parcel-level `A`/`M` metrics then fits
genpca on stacked runs. See voxel-level
[`fit_subject_genpca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca.md)
for an overview.

## Usage

``` r
fit_subject_genpca_parcel(
  cnv_list,
  parc_vol = NULL,
  parcel_coords = NULL,
  adjacency = NULL,
  gm_vol = NULL,
  wm_vol = NULL,
  csf_vol = NULL,
  gm_p = NULL,
  wm_p = NULL,
  csf_p = NULL,
  wm_parcels = NULL,
  FD_list = NULL,
  DVARS_list = NULL,
  k = 20,
  p_ar = 1L,
  lambda_s = 0.5,
  lambda_t = 0.3,
  knn_k = 6,
  knn_sigma = 2.5
)
```

## Arguments

- cnv_list:

  List of `ClusteredNeuroVec` runs.

- parc_vol, parcel_coords, adjacency:

  Ways to define the parcel graph: provide exactly one.

- gm_vol, wm_vol, csf_vol:

  Optional tissue maps used if `gm_p/wm_p/csf_p` are not provided.

- gm_p, wm_p, csf_p:

  Optional parcel-level tissue summaries (avoids recomputation).

- wm_parcels:

  Optional WM parcel indices for AR whitener estimation.

- FD_list, DVARS_list:

  Optional lists of per-run motion metrics. Each list should contain
  numeric vectors (one per run in cnv_list). FD is framewise
  displacement (mm), DVARS is temporal derivative of signal variance.
  Used to down-weight high-motion frames. See
  [`make_frame_weights`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md)
  for detailed definitions of these metrics.

- k:

  Number of components to extract (default 20).

- p_ar:

  AR order for temporal whitening (default 1).

- lambda_s:

  Spatial regularization weight (default 0.5).

- lambda_t:

  Temporal regularization weight (default 0.3).

- knn_k:

  Number of nearest neighbors for parcel graph (default 6).

- knn_sigma:

  Gaussian kernel bandwidth for parcel graph (default 2.5).

## Value

A `genpca` object with parcel-level loadings.

## Details

This is the parcel-level analogue of
[`fit_subject_genpca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca.md),
operating on `ClusteredNeuroVec` data where voxels have been aggregated
into parcels. The workflow is:

1.  Build parcel graph Laplacian from adjacency, coordinates, or labels

2.  Compute parcel tSNR across runs

3.  Aggregate tissue maps to parcels (if not provided)

4.  Construct column metric A with spatial and tissue weights

5.  Build per-run row metrics M with AR whitening and motion weights

6.  Stack runs temporally and fit generalized PCA

Parcellation reduces dimensionality and computation while preserving
spatial structure through the parcel graph.

## See also

[`fit_subject_genpca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca.md)
for voxel-level analysis

Other main fitting:
[`fit_subject_genpca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca.md),
[`fit_subject_metapca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_metapca.md)

## Examples

``` r
# \donttest{
# Fit genpca on parcellated data with tissue weights
fit <- fit_subject_genpca_parcel(cnv_list, parc_vol = parcellation,
                                 gm_vol = gm, wm_vol = wm, csf_vol = csf,
                                 k = 20)
#> Error: object 'cnv_list' not found
# }
```
