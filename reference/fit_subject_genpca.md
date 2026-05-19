# Fit generalized PCA across multiple runs (voxel level)

Builds a voxel-level column metric `A` and per-run row metrics `M_r`,
stacks runs in time, and fits `genpca` with `A` and block-diagonal `M`.
Adds `scores` (projected time scores) for convenience.

## Usage

``` r
fit_subject_genpca(
  nv_list,
  mask_vol,
  gm_vol,
  wm_vol,
  csf_vol,
  FD_list = NULL,
  DVARS_list = NULL,
  k = 20,
  p_ar = 1L,
  lambda_s = 0.5,
  lambda_t = 0.3,
  lap_k = 6,
  lap_sigma = 2.5
)
```

## Arguments

- nv_list:

  List of `NeuroVec` runs.

- mask_vol:

  `NeuroVol` binary mask.

- gm_vol, wm_vol, csf_vol:

  `NeuroVol` tissue probability maps (same space as mask).

- FD_list:

  Optional list of per-run framewise displacement (FD) values. Should
  contain numeric vectors, one per run in nv_list. FD measures head
  motion in mm between consecutive frames. Used to down-weight
  high-motion frames. See
  [`make_frame_weights`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md)
  for detailed definition.

- DVARS_list:

  Optional list of per-run DVARS (temporal derivative of signal
  variance) values. Should contain numeric vectors, one per run in
  nv_list. DVARS measures frame-to-frame signal change. Used to
  down-weight frames with excessive signal variation. See
  [`make_frame_weights`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md)
  for detailed definition.

- k:

  Number of components to fit (default 20). Start with modest values
  (5-10) for exploration, increase based on rank selection diagnostics.

- p_ar:

  AR order for row whitener (default 1). Use 0 to skip AR modeling,
  higher values (2-3) for data with stronger temporal autocorrelation.

- lambda_s:

  Spatial regularization weight (default 0.5). Controls smoothness via
  the graph Laplacian; higher values enforce more spatial coherence.

- lambda_t:

  Temporal regularization weight (default 0.3). Controls temporal
  smoothness; higher values penalize rapid fluctuations.

- lap_k:

  Number of nearest neighbors for Laplacian construction (default 6).
  Typically 6-8 for 3D connectivity; higher values increase spatial
  coupling.

- lap_sigma:

  Gaussian kernel bandwidth for edge weights (default 2.5). Controls the
  decay of spatial influence with distance.

## Value

A `genpca` object with loadings `v` and added `scores` (stacked time x
k).

## Details

Validates that all `NeuroVec` runs share the same spatial `NeuroSpace`
as the mask and tissue maps. Constructs a Laplacian `L` via
`make_laplacian`, a tissue/tSNR-weighted `A` via `build_spatial_metric`,
and per-run `M_r` via `make_M_for_run` using optional FD/DVARS. Rows
(time) are stacked across runs and
[`genpca::genpca`](https://bbuchsbaum.github.io/genpca/reference/genpca.html)
is called with `ncomp = k`. The result gains a `scores` component via
[`multivarious::scores`](https://bbuchsbaum.github.io/multivarious/reference/scores.html).

## See also

[`fit_subject_genpca_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca_parcel.md)
for parcel-level analysis,
[`fit_subject_metapca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_metapca.md)
for combining multiple fits

Other main fitting:
[`fit_subject_genpca_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca_parcel.md),
[`fit_subject_metapca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_metapca.md)

## Examples

``` r
# \donttest{
fit <- fit_subject_genpca(nv_list, mask, gm, wm, csf, k = 5)
#> Error: object 'nv_list' not found
# }
```
