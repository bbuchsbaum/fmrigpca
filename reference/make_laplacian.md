# Build a spatial graph Laplacian from a NeuroVol mask

Constructs a sparse, symmetric voxel–wise graph Laplacian for a masked
volume using k–nearest neighbors with a Gaussian kernel. Optionally
returns the normalized Laplacian and scales by the largest eigenvalue
for numerical stability.

## Usage

``` r
make_laplacian(mask_vol, k = 6, sigma = 2.5, normalized = TRUE)
```

## Arguments

- mask_vol:

  A
  [`neuroim2::NeuroVol`](https://bbuchsbaum.github.io/neuroim2/reference/NeuroVol.html)
  logical/numeric mask defining in-mask voxels.

- k:

  Integer; number of neighbors in the k–NN graph (default 6). Common
  choices: 6 for face connectivity, 18 for face+edge, 26 for full 3D.

- sigma:

  Numeric; Gaussian kernel bandwidth for edge weights (default 2.5).
  Smaller values create sparser graphs; larger values smooth more
  broadly.

- normalized:

  Logical; return the normalized Laplacian if `TRUE` (default).
  Normalized version has better numerical properties and bounded
  spectrum.

## Value

A sparse, symmetric `Matrix` of size `n_voxels x n_voxels` representing
the (optionally normalized) graph Laplacian over in-mask voxels.

## Details

Voxels inside the mask are converted to integer coordinates (i,j,k),
then a k–NN graph is built with edge weights `exp(-d^2 / (2*sigma^2))`
using the `adjoin` package (`graph_weights`/`adjacency`) when available,
falling back to an internal implementation otherwise. The unnormalized
Laplacian is `L = D - W`. The normalized variant is
`L = I - D^{-1/2} W D^{-1/2}`. In both cases, `L` is divided by its
largest eigenvalue when positive and finite, keeping the spectrum within
a stable range (e.g., 0 to 2 for the normalized Laplacian).

Edge cases are handled: a single-voxel mask yields a 1x1 identity; empty
masks are rejected. If `k` is too large, it is clipped to `n-1`
neighbors.

## See also

[`make_parcel_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_parcel_laplacian.md)
for parcel-level Laplacian construction

Other spatial metrics:
[`build_spatial_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric.md),
[`build_spatial_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric_parcel.md),
[`compute_tsnr()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr.md),
[`compute_tsnr_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr_parcel.md),
[`make_parcel_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_parcel_laplacian.md)

## Examples

``` r
# \donttest{
dims <- c(6,6,4)
arr  <- array(FALSE, dims); arr[3:4,3:4,2:3] <- TRUE
mv   <- neuroim2::NeuroVol(arr, neuroim2::NeuroSpace(dims, c(2,2,2)))
L    <- make_laplacian(mv, k = 4, sigma = 2, normalized = TRUE)
Matrix::isSymmetric(L)
#> [1] TRUE
# }
```
