# Temporal SNR (tSNR) per voxel across runs

Computes per-voxel temporal signal-to-noise ratio across a list of
`NeuroVec` runs by aggregating first and second moments over time and
runs.

## Usage

``` r
compute_tsnr(nv_list)
```

## Arguments

- nv_list:

  List of
  [`neuroim2::NeuroVec`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVec-class.html)
  runs.

## Value

Numeric vector of length V (voxels) with tSNR values.

## Details

For each voxel, the mean and variance are computed over concatenated
time points across runs. tSNR is `max(1e-6, mean) / max(1e-6, sd)`, then
lower- bounded at `1e-6` to avoid degenerate weights.

## See also

[`compute_tsnr_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr_parcel.md)
for parcel-level tSNR computation

Other spatial metrics:
[`build_spatial_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric.md),
[`build_spatial_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric_parcel.md),
[`compute_tsnr_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr_parcel.md),
[`make_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_laplacian.md),
[`make_parcel_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_parcel_laplacian.md)

## Examples

``` r
# \donttest{
# Create simple test data
dims <- c(10, 10, 5, 20)  # 4D: x, y, z, time
space <- neuroim2::NeuroSpace(dims,
                               spacing = c(1, 1, 1),
                               origin = c(0, 0, 0))
mask <- array(TRUE, dims[1:3])  # 3D mask

# Create two runs with random data
data1 <- array(rnorm(prod(dims)), dims)
data2 <- array(rnorm(prod(dims)), dims)
nv1 <- neuroim2::NeuroVec(data1, space, mask = mask)
nv2 <- neuroim2::NeuroVec(data2, space, mask = mask)

tsnr <- compute_tsnr(list(nv1, nv2))
# }
```
