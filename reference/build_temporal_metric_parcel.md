# Build temporal covariance metric for parcellated fMRI data

Constructs a temporal covariance metric (row metric M) for parcellated
fMRI data in `ClusteredNeuroVec` format. This is the parcel-level
analogue of
[`build_temporal_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric.md)
for voxel data.

## Usage

``` r
build_temporal_metric_parcel(
  cnv_run,
  wm_parcels = NULL,
  p = 1L,
  FD = NULL,
  DVARS = NULL,
  lambda_t = 0.3,
  ridge = 1e-06
)
```

## Arguments

- cnv_run:

  A `ClusteredNeuroVec` run.

- wm_parcels:

  Optional integer indices of WM parcels for AR estimation.

- p:

  AR order (default 1). Higher orders model longer-range dependencies.

- FD, DVARS:

  Optional motion metrics for the run. FD is framewise displacement
  (mm), DVARS is temporal derivative of signal variance. See
  [`make_frame_weights`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md)
  for detailed definitions.

- lambda_t:

  Temporal penalty weight (default 0.3). Higher values enforce more
  temporal smoothness.

- ridge:

  Small ridge added to ensure positive definiteness (default 1e-6).

## Value

Symmetric positive definite `Matrix` (T x T).

## Details

The metric is `M = W^{1/2} Q' H Q W^{1/2} + ridge * I`, where:

- `Q` whitens the temporal covariance from AR(p) estimation

- `H` is the second-difference temporal smoothness penalty

- `W` is diagonal frame weighting from motion parameters

This parcel-level implementation uses parcel-averaged time series for AR
estimation, providing computational efficiency for high-dimensional
data.

## See also

[`build_temporal_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric.md)
for voxel-level temporal metric construction

Other temporal metrics:
[`build_temporal_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric.md),
[`estimate_ar_whitener()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener.md),
[`estimate_ar_whitener_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener_parcel.md),
[`make_frame_weights()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md),
[`make_temporal_penalty()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_temporal_penalty.md)

## Examples

``` r
# \donttest{
# Build temporal metric with AR(1) whitening and motion regressors
M <- build_temporal_metric_parcel(cnv_run, wm_parcels = c(1,2,3),
                                  FD = fd_vals, DVARS = dvars_vals)
#> Error: object 'cnv_run' not found
# }
```
