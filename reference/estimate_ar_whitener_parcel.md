# Estimate AR(p) whitener from parcel-level data

Parcel analogue of
[`estimate_ar_whitener()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener.md)
that uses the parcel mean or optional `wm_parcels` subset to estimate
AR(p) via Yule–Walker.

## Usage

``` r
estimate_ar_whitener_parcel(cnv_run, wm_parcels = NULL, p = 1L)
```

## Arguments

- cnv_run:

  A `ClusteredNeuroVec` run.

- wm_parcels:

  Optional integer indices of WM parcels for more robust AR estimation.

- p:

  Integer AR order (default 1). Use higher orders for data with stronger
  temporal autocorrelation.

## Value

List with components:

- Q:

  Inverse Cholesky factor (whitening matrix)

- Sigma:

  Regularized AR covariance matrix

- phi:

  AR coefficients of length p

## Details

The time series is formed as the mean across parcels, optionally
restricted to white-matter parcels if `wm_parcels` indices are provided.
The series is standardized and AR(p) coefficients are estimated via
[`stats::ar.yw`](https://rdrr.io/r/stats/ar.html). The implied Toeplitz
covariance matrix is regularized to ensure positive definiteness.
Returns the whitener `Q = chol(Sigma)^{-1}` for use in the row metric
construction.

## See also

[`estimate_ar_whitener()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener.md)
for voxel-level AR whitening

Other temporal metrics:
[`build_temporal_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric.md),
[`build_temporal_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric_parcel.md),
[`estimate_ar_whitener()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener.md),
[`make_frame_weights()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md),
[`make_temporal_penalty()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_temporal_penalty.md)

## Examples

``` r
# \donttest{
# Estimate AR(1) whitener using WM parcels
arw <- estimate_ar_whitener_parcel(cnv_run, wm_parcels = c(1,2,3), p = 1)
#> Error: object 'cnv_run' not found
# }
```
