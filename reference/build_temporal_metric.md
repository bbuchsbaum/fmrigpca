# Build temporal covariance metric for generalized PCA

Constructs a temporal covariance metric (row metric M) that encodes
temporal dependencies, smoothness constraints, and data quality weights
for a single fMRI run. This metric is used as the row-wise (temporal)
metric in generalized PCA to properly account for temporal structure in
fMRI data.

## Usage

``` r
build_temporal_metric(
  nv_run,
  wm_mask = NULL,
  p = 1L,
  FD = NULL,
  DVARS = NULL,
  lambda_t = 0.3,
  ridge = 1e-06
)
```

## Arguments

- nv_run:

  A
  [`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
  (single run) containing fMRI time series.

- wm_mask:

  Optional `NeuroVol` white matter mask. When provided, AR parameters
  are estimated from the mean white matter signal, which typically has
  less neural signal and better represents noise autocorrelation.

- p:

  Integer AR order (default 1). Common choices:

  - 0: No AR modeling (independent time points)

  - 1: First-order autoregression (standard, captures ~70% of
    autocorrelation)

  - 2-3: Higher orders for data with stronger temporal dependencies

- FD:

  Optional numeric vector of framewise displacement values (mm). Length
  must equal the number of time points. Typically computed during motion
  correction preprocessing.

- DVARS:

  Optional numeric vector of DVARS values (signal derivative RMS).
  Length must equal the number of time points. Typically computed from
  preprocessed BOLD signal.

- lambda_t:

  Temporal smoothness penalty weight (default 0.3). Range 0-1:

  - 0: No temporal smoothing

  - 0.3: Standard smoothing (default)

  - 0.5-0.7: Strong smoothing for high-motion or noisy data

  - 1.0: Maximum smoothing

- ridge:

  Small positive value added to diagonal for numerical stability
  (default 1e-6). Increase if encountering singular matrix errors.

## Value

A symmetric positive definite `Matrix` of size `T x T` encoding temporal
covariance structure for use in generalized PCA.

## Details

### Purpose in Generalized PCA

In standard PCA, all time points are treated as independent and equally
reliable. This is inappropriate for fMRI data where:

- Consecutive time points are correlated (temporal autocorrelation)

- Some time points are corrupted by motion artifacts

- Neural signals vary smoothly over time

This function builds a metric M that encodes these properties, allowing
generalized PCA to find components that respect the temporal structure
of fMRI data.

### Mathematical Formulation

The temporal metric combines three key components:

`M = W^{1/2} Q' H Q W^{1/2} + ridge * I`

where:

- `Q`: AR(p) whitening matrix that removes temporal autocorrelation

- `H`: Temporal smoothness penalty (second-difference operator)

- `W`: Diagonal weights based on motion quality (FD/DVARS)

- `ridge`: Small regularization for numerical stability

### Components Explained

**AR Whitening (Q)**: fMRI time series exhibit strong autocorrelation
(typically 0.3-0.5 at lag 1). The AR whitener transforms the data to
make time points approximately independent, preventing
overrepresentation of slow temporal drifts in the components.

**Temporal Smoothness (H)**: Neural signals change smoothly over time.
The penalty `H = I + λ_t * D2' D2` (where D2 is the second-difference
operator) encourages components with smooth temporal profiles while
allowing for dynamic changes.

**Motion Weighting (W)**: Time points with high motion (FD) or global
signal changes (DVARS) are down-weighted rather than discarded. This
soft censoring approach maintains temporal continuity while reducing
artifact influence. See
[`make_frame_weights()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md)
for details.

### Usage Scenarios

- **Standard**: `p=1, lambda_t=0.3` with both FD and DVARS

- **High motion data**: Increase lambda_t (0.5-0.7) for more smoothing

- **Task fMRI**: Consider p=0 if modeling task regressors separately

- **Resting state**: Always use p≥1 to account for autocorrelation

## See also

- [`build_temporal_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric_parcel.md)
  for parcel-level temporal metrics

- [`estimate_ar_whitener()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener.md)
  for AR parameter estimation details

- [`make_temporal_penalty()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_temporal_penalty.md)
  for smoothness penalty construction

- [`make_frame_weights()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md)
  for motion-based weighting

Other temporal metrics:
[`build_temporal_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric_parcel.md),
[`estimate_ar_whitener()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener.md),
[`estimate_ar_whitener_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener_parcel.md),
[`make_frame_weights()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_frame_weights.md),
[`make_temporal_penalty()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_temporal_penalty.md)

## Examples

``` r
# \donttest{
# Basic usage with AR(1) whitening
time_points <- 12
dims <- c(4, 4, 3, time_points)  # 4D: x, y, z, time
space <- neuroim2::NeuroSpace(dims,
                               spacing = c(1, 1, 1),
                               origin = c(0, 0, 0))
mask <- array(TRUE, dims[1:3])  # 3D mask
data <- array(rnorm(prod(dims)), dims)
nv_run <- neuroim2::NeuroVec(data, space, mask = mask)
wm_vol <- neuroim2::NeuroVol(array(runif(prod(dims[1:3])), dims[1:3]),
                             neuroim2::NeuroSpace(dims[1:3], c(1, 1, 1)))
fd_values <- runif(time_points, 0, 0.4)
dvars_values <- runif(time_points, 0, 1)
M <- build_temporal_metric(nv_run, p = 1L)

# With motion correction
M <- build_temporal_metric(nv_run, wm_mask = wm_vol, 
                           FD = fd_values, DVARS = dvars_values,
                           lambda_t = 0.5)

# No AR modeling for task fMRI with explicit task regressors
M <- build_temporal_metric(nv_run, p = 0, lambda_t = 0.3)
# }
```
