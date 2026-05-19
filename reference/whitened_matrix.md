# Form a whitened matrix for genpca geometry

Applies square-root scalings of the row and column metrics to obtain the
whitened matrix `W = R_M^T X R_A`, where `A = R_A^T R_A` and
`M = R_M^T R_M`.

## Usage

``` r
whitened_matrix(X, A, M)
```

## Arguments

- X:

  Numeric matrix of size T x V (time x voxels/parcels).

- A:

  Column metric (V x V), symmetric positive definite.

- M:

  Row metric (T x T), symmetric positive definite.

## Value

The whitened matrix `W` as a `Matrix`.

## See also

[`whitened_svd()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_svd.md)
for computing SVD of the whitened matrix

Other subspace analysis:
[`procrustes_distance()`](https://bbuchsbaum.github.io/fmrigpca/reference/procrustes_distance.md),
[`subspace_principal_angles()`](https://bbuchsbaum.github.io/fmrigpca/reference/subspace_principal_angles.md),
[`whitened_svd()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_svd.md)

## Examples

``` r
# \donttest{
set.seed(1)
X <- matrix(rnorm(200), nrow = 20, ncol = 10)
A <- Matrix::Diagonal(ncol(X))
M <- Matrix::Diagonal(nrow(X))
W <- whitened_matrix(X, A, M)
# }
```
