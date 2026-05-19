# Build spatial column metric A from tissue maps, Laplacian, and tSNR

Constructs a positive definite spatial column metric `A` that encodes
spatial smoothness (graph Laplacian) and voxelwise weights derived from
tissue probabilities and temporal SNR. Adds a small ridge for numerical
stability.

## Usage

``` r
build_spatial_metric(
  gm_vol,
  wm_vol,
  csf_vol,
  L,
  tsnr = NULL,
  mask_idx = NULL,
  alpha = 1,
  beta = 0.5,
  gamma = 1,
  lambda_s = 0.5,
  tau = 1e-06
)
```

## Arguments

- gm_vol, wm_vol, csf_vol:

  `NeuroVol` tissue probability maps.

- L:

  Laplacian `Matrix` over in-mask voxels.

- tsnr:

  Optional numeric tSNR vector over in-mask voxels.

- mask_idx:

  Integer indices of in-mask voxels (length equals `nrow(L)`).

- alpha, beta, gamma:

  Exponents for GM, tSNR, and (WM+CSF) weighting.

- lambda_s:

  Spatial regularization weight multiplying `L`.

- tau:

  Small ridge added to the diagonal (default 1e-6).

## Value

Symmetric positive definite `Matrix` of size `V_in x V_in`.

## Details

Let `w = gm^alpha * tsnr^beta * (wm + csf)^{-gamma}` on in-mask voxels.
Define `S = I + lambda_s * L` where `L` is a (normalized) Laplacian.
Then `A = diag(sqrt(w)) S diag(sqrt(w)) + tau * I`. The `mask_idx`
selects in-mask voxel indices and must align with `nrow(L)`. If `tsnr`
is `NULL`, it defaults to 1. Tissue maps must share space with the mask.

## See also

[`build_spatial_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric_parcel.md)
for parcel-level column metric construction

Other spatial metrics:
[`build_spatial_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric_parcel.md),
[`compute_tsnr()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr.md),
[`compute_tsnr_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr_parcel.md),
[`make_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_laplacian.md),
[`make_parcel_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_parcel_laplacian.md)

## Examples

``` r
# \donttest{
# Create small example volumes
dims <- c(5, 5, 3)
space <- neuroim2::NeuroSpace(dims, c(1, 1, 1))

# Create tissue probability maps
gm_data <- array(runif(prod(dims), 0.5, 1), dims)
wm_data <- array(runif(prod(dims), 0, 0.3), dims)
csf_data <- array(runif(prod(dims), 0, 0.2), dims)

gm <- neuroim2::NeuroVol(gm_data, space)
wm <- neuroim2::NeuroVol(wm_data, space)
csf <- neuroim2::NeuroVol(csf_data, space)

# Create mask and Laplacian
mask <- neuroim2::NeuroVol(gm_data > 0.3, space)
L <- make_laplacian(mask, k = 6)
mask_idx <- which(neuroim2::values(mask) > 0)

A <- build_spatial_metric(gm, wm, csf, L, mask_idx = mask_idx)
# }
```
