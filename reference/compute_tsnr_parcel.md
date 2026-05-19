# Parcel-level tSNR across runs (ClusteredNeuroVec)

Computes parcel-wise temporal signal-to-noise ratio across
`ClusteredNeuroVec` runs by aggregating means and variances over time.

## Usage

``` r
compute_tsnr_parcel(cnv_list)
```

## Arguments

- cnv_list:

  List of
  [`neuroim2::ClusteredNeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/ClusteredNeuroVec.html)
  runs.

## Value

Numeric vector of length P with parcel tSNR values.

## Details

For each parcel, the mean and variance are computed over concatenated
time points across all runs. tSNR is calculated as `mean / sd`, with
both numerator and denominator lower-bounded at `1e-6` to avoid division
by zero or degenerate weights. This is the parcel-level analogue of
[`compute_tsnr()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr.md)
for voxel data.

## See also

[`compute_tsnr()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr.md)
for voxel-level tSNR computation

Other spatial metrics:
[`build_spatial_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric.md),
[`build_spatial_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric_parcel.md),
[`compute_tsnr()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr.md),
[`make_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_laplacian.md),
[`make_parcel_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_parcel_laplacian.md)

## Examples

``` r
# \donttest{
# Compute parcel tSNR across multiple runs
tsnr_p <- compute_tsnr_parcel(cnv_list)
#> Error: object 'cnv_list' not found
# }
```
