# Truncated SVD of the whitened matrix

Computes singular values (optionally truncated) of the whitened matrix.
Uses `RSpectra` when available for efficiency in the truncated case.

## Usage

``` r
whitened_svd(X, A, M, k = NULL)
```

## Arguments

- X, A, M:

  As in
  [`whitened_matrix()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_matrix.md).

- k:

  Optional integer number of singular values to compute.

## Value

A list with `d` (singular values) and `u`,`v` set to `NULL`.

## See also

[`whitened_matrix()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_matrix.md)
for the underlying whitening operation

Other subspace analysis:
[`procrustes_distance()`](https://bbuchsbaum.github.io/fmrigpca/reference/procrustes_distance.md),
[`subspace_principal_angles()`](https://bbuchsbaum.github.io/fmrigpca/reference/subspace_principal_angles.md),
[`whitened_matrix()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_matrix.md)

## Examples

``` r
# \donttest{
set.seed(1)
X <- matrix(rnorm(300), nrow = 30, ncol = 10)
A <- Matrix::Diagonal(ncol(X))
M <- Matrix::Diagonal(nrow(X))
ws <- whitened_svd(X, A, M, k = 5)
# }
```
