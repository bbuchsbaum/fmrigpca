# Build parcel-level spatial column metric A from aggregated tissues and Laplacian

Forms a parcel-level spatial metric `A` analogous to voxel-level
`build_spatial_metric`, using aggregated tissue probabilities and
optional parcel tSNR.

## Usage

``` r
build_spatial_metric_parcel(
  gm_p,
  wm_p,
  csf_p,
  Lp,
  tsnr_p = NULL,
  alpha = 1,
  beta = 0.5,
  gamma = 1,
  lambda_s = 0.5,
  tau = 1e-06
)
```

## Arguments

- gm_p, wm_p, csf_p:

  Numeric vectors of length P with parcel tissue values (typically means
  of voxel-wise tissue probabilities within each parcel).

- Lp:

  Parcel Laplacian (P x P) from
  [`make_parcel_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_parcel_laplacian.md).

- tsnr_p:

  Optional parcel tSNR vector from
  [`compute_tsnr_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr_parcel.md).

- alpha:

  Exponent for GM weighting (default 1.0). Higher values give more
  weight to gray matter parcels.

- beta:

  Exponent for tSNR weighting (default 0.5). Controls influence of
  temporal signal quality.

- gamma:

  Exponent for (WM+CSF) down-weighting (default 1.0). Higher values
  penalize non-gray matter more strongly.

- lambda_s:

  Spatial regularization weight (default 0.5). Controls the strength of
  spatial smoothing via the Laplacian.

- tau:

  Small ridge on the diagonal (default 1e-6) for numerical stability.

## Value

Symmetric positive definite `Matrix` (P x P).

## Details

Constructs weights `w = gm^alpha * tsnr^beta * (wm + csf)^{-gamma}` for
each parcel, then forms
`A = diag(sqrt(w)) * (I + lambda_s * Lp) * diag(sqrt(w)) + tau * I`.
This encodes both spatial smoothness (via the Laplacian) and parcel-wise
quality weights derived from tissue composition and temporal SNR. The
result is a positive definite metric suitable for generalized PCA.

## See also

[`build_spatial_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric.md)
for voxel-level column metric construction

Other spatial metrics:
[`build_spatial_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric.md),
[`compute_tsnr()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr.md),
[`compute_tsnr_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr_parcel.md),
[`make_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_laplacian.md),
[`make_parcel_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_parcel_laplacian.md)

## Examples

``` r
# \donttest{
# Build parcel column metric with tissue weights and spatial smoothing
gm_p <- c(0.8, 0.6, 0.2)
wm_p <- c(0.15, 0.25, 0.6)
csf_p <- c(0.05, 0.15, 0.2)
tsnr_vals <- c(50, 40, 30)
adj <- Matrix::Matrix(matrix(c(0, 1, 0,
                               1, 0, 1,
                               0, 1, 0), nrow = 3, byrow = TRUE))
deg <- Matrix::Diagonal(x = c(1, 2, 1))
Lp <- deg - adj
A <- build_spatial_metric_parcel(gm_p, wm_p, csf_p, Lp, tsnr_p = tsnr_vals,
                                 lambda_s = 0.5)
# }
```
