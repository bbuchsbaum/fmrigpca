# Frame-wise robust weights from FD and DVARS

Converts motion summaries (FD) and DVARS into diagonal frame weights in
(0,1\], down-weighting high-motion or high-DVARS frames. If both inputs
are `NULL`, returns a `1x1` identity (caller may expand to `T x T`).

## Usage

``` r
make_frame_weights(FD = NULL, DVARS = NULL, fd_thresh = 0.5, dvars_z = 3)
```

## Arguments

- FD:

  Optional numeric vector of framewise displacement values in mm.
  Typically computed from motion correction parameters during
  preprocessing.

- DVARS:

  Optional numeric vector of DVARS values (signal RMS derivative).
  Typically computed from preprocessed BOLD time series.

- fd_thresh:

  Threshold above which FD begins to down-weight (default 0.5mm). Common
  choices: 0.2mm (strict), 0.5mm (standard), 1.0mm (lenient).

- dvars_z:

  Z-score threshold for DVARS down-weighting (default 3.0). Frames with
  DVARS z-scores above this are down-weighted; 2.5-3.5 typical.

## Value

A diagonal `Matrix` of weights.

## Details

### Framewise Displacement (FD)

FD quantifies head motion between consecutive fMRI volumes as defined by
Power et al. (2012). It is calculated as the sum of absolute values of
the temporal derivatives of the six rigid-body motion parameters (3
translations and 3 rotations):

`FD = |Δx| + |Δy| + |Δz| + |Δα| + |Δβ| + |Δγ|`

where rotational parameters (α, β, γ) are converted to millimeters by
calculating arc length displacement on a sphere with radius 50mm. This
conversion ensures all parameters are on the same scale. FD values
represent the total displacement in mm from one volume to the next.

### DVARS (Derivative of VARS)

DVARS measures the rate of change of BOLD signal intensity across the
whole brain between consecutive time points (Smyser et al., 2010; Power
et al., 2012). It is calculated as the root mean square of the temporal
derivative:

`DVARS_t = sqrt(mean((x_t - x_{t-1})^2))`

where x represents voxel intensities. DVARS is sensitive to both motion
artifacts and physiological noise. High DVARS values indicate sudden
global signal changes that often correlate with head motion.

### Weight Calculation

This function converts these metrics to frame weights:

- FD weights: `w_FD = 1 / (1 + max(0, FD - fd_thresh))`

- DVARS weights: `w_DVARS = 1 / (1 + max(0, z - dvars_z))` where z is
  the z-score

- Final weight: `w = max(0.001, w_FD * w_DVARS)`

The resulting diagonal weight matrix down-weights problematic frames
while preserving data continuity, as complete censoring can introduce
discontinuities in the time series.

## References

Power JD, Barnes KA, Snyder AZ, Schlaggar BL, Petersen SE (2012).
"Spurious but systematic correlations in functional connectivity MRI
networks arise from subject motion." NeuroImage, 59(3), 2142-2154.

Smyser CD, Inder TE, Shimony JS, Hill JE, Degnan AJ, Snyder AZ, Neil JJ
(2010). "Longitudinal analysis of neural network development in preterm
infants." Cerebral Cortex, 20(12), 2852-2862.

## See also

Other temporal metrics:
[`build_temporal_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric.md),
[`build_temporal_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric_parcel.md),
[`estimate_ar_whitener()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener.md),
[`estimate_ar_whitener_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener_parcel.md),
[`make_temporal_penalty()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_temporal_penalty.md)

## Examples

``` r
# Create frame weights from motion parameters
W <- make_frame_weights(FD = runif(30, 0, 0.5), DVARS = runif(30, 0.5, 1.5))
Matrix::diag(W)[1:5]  # First 5 frame weights
#> [1] 1 1 1 1 1
```
