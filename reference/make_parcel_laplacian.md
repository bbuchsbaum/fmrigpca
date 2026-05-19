# Build a parcel-level Laplacian

Constructs a parcel adjacency-based Laplacian either from a supplied
adjacency matrix, parcel centroids (k–NN), or a label volume.
Normalization and spectral scaling mirror the voxel-level Laplacian.

## Usage

``` r
make_parcel_laplacian(
  adjacency = NULL,
  parcel_coords = NULL,
  parc_vol = NULL,
  k = 6,
  sigma = 2.5,
  normalized = TRUE
)
```

## Arguments

- adjacency:

  Optional symmetric adjacency/weight matrix (P x P).

- parcel_coords:

  Optional P x 3 matrix of parcel centroids.

- parc_vol:

  Optional `NeuroVol` of parcel labels (used to compute centroids).

- k, sigma:

  KNN neighbors and Gaussian bandwidth.

- normalized:

  Logical; return normalized Laplacian.

## Value

Sparse symmetric `Matrix` (P x P) Laplacian.

## See also

[`make_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_laplacian.md)
for voxel-level Laplacian construction

Other spatial metrics:
[`build_spatial_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric.md),
[`build_spatial_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric_parcel.md),
[`compute_tsnr()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr.md),
[`compute_tsnr_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/compute_tsnr_parcel.md),
[`make_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_laplacian.md)

## Examples

``` r
# \donttest{
Lp <- make_parcel_laplacian(parc_vol = parc, k = 6, sigma = 2)
#> Error: object 'parc' not found
# }
```
