# Principal angles between two subspaces

Computes principal angles between the column spaces spanned by `U` and
`V`.

## Usage

``` r
subspace_principal_angles(U, V)
```

## Arguments

- U, V:

  Matrices with the same number of rows.

## Value

Numeric vector of angles in radians.

## See also

[`procrustes_distance()`](https://bbuchsbaum.github.io/fmrigpca/reference/procrustes_distance.md)
for an alternative subspace distance metric

Other subspace analysis:
[`procrustes_distance()`](https://bbuchsbaum.github.io/fmrigpca/reference/procrustes_distance.md),
[`whitened_matrix()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_matrix.md),
[`whitened_svd()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_svd.md)

## Examples

``` r
# \donttest{
set.seed(1)
U <- matrix(rnorm(50), nrow = 10)
V <- matrix(rnorm(50), nrow = 10)
theta <- subspace_principal_angles(U, V)
# }
```
