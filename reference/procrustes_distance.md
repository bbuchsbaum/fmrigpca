# Procrustes distance between two subspaces

Computes the orthogonal Procrustes distance between subspaces spanned by
`U` and `V` using singular values of `U'V`.

## Usage

``` r
procrustes_distance(U, V)
```

## Arguments

- U, V:

  Basis matrices for the two subspaces.

## Value

Nonnegative scalar distance.

## See also

[`subspace_principal_angles()`](https://bbuchsbaum.github.io/fmrigpca/reference/subspace_principal_angles.md)
for computing angles between subspaces

Other subspace analysis:
[`subspace_principal_angles()`](https://bbuchsbaum.github.io/fmrigpca/reference/subspace_principal_angles.md),
[`whitened_matrix()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_matrix.md),
[`whitened_svd()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_svd.md)

## Examples

``` r
# \donttest{
set.seed(1)
U <- matrix(rnorm(50), nrow = 10)
V <- matrix(rnorm(50), nrow = 10)
d <- procrustes_distance(U, V)
# }
```
