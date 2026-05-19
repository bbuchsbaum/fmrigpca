# Estimate a design-free AR(p) whitener from a run

Estimates an autoregressive whitener for a single `NeuroVec` run using
the global mean (optionally restricted to WM) and Yule–Walker AR
fitting. Returns the inverse Cholesky `Q`, the covariance `Sigma`, and
AR coefficients `phi`.

## Usage

``` r
estimate_ar_whitener(nv_run, wm_mask = NULL, p = 1L)
```

## Arguments

- nv_run:

  A
  [`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
  (single run) with dimensions V x T.

- wm_mask:

  Optional `NeuroVol` to restrict the mean to WM voxels.

- p:

  Integer AR order (default 1). Use 0 to skip AR modeling.

## Value

A list with components `Q` (inverse Cholesky), `Sigma` (covariance), and
`phi` (AR coefficients of length `p`).

## Details

The time series is formed as the mean across voxels, optionally within a
provided white-matter mask. It is standardized, and AR(p) coefficients
are estimated via [`stats::ar.yw`](https://rdrr.io/r/stats/ar.html). The
implied Toeplitz covariance `Sigma` is regularized to positive
definiteness (via a ridge and
[`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html) if
needed). The whitener `Q = chol(Sigma)^{-1}` is returned for use in the
row metric. If `p = 0`, an identity whitener is returned.

## See also

[`estimate_ar_whitener_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener_parcel.md)
for parcel-level AR whitening

Other temporal metrics:
[`build_temporal_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric.md),
[`build_temporal_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric_parcel.md),
[`estimate_ar_whitener_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener_parcel.md),
[`make_frame_weights()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md),
[`make_temporal_penalty()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_temporal_penalty.md)

## Examples

``` r
# \donttest{
nv_run <- neuroim2::simulate_fmri(mask, n_time = 50, seed = 1)
#> Error: object 'mask' not found
arw   <- estimate_ar_whitener(nv_run, p = 1L)
#> Error: object 'nv_run' not found
# }
```
