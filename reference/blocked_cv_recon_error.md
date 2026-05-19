# Blocked cross-validated reconstruction error under genpca geometry

Evaluates reconstruction error by holding out contiguous time blocks and
projecting hold-out data onto the `k`-dimensional column space learned
on the training block, all under the whitened geometry.

## Usage

``` r
blocked_cv_recon_error(X, A, M, k, nfold = 5L, block_frac = 0.2)
```

## Arguments

- X, A, M:

  Data and metrics as in
  [`whitened_matrix()`](https://bbuchsbaum.github.io/fmrigpca/reference/whitened_matrix.md).

- k:

  Number of components.

- nfold:

  Number of folds (hold-out blocks).

- block_frac:

  Fraction of time points per hold-out block.

## Value

Data frame with per-fold errors and the mean.

## See also

[`choose_rank_gd()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_gd.md)
for Gavish-Donoho approach,
[`choose_rank_pa()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_pa.md)
for parallel analysis approach

Other rank selection:
[`choose_rank_gd()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_gd.md),
[`choose_rank_pa()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_pa.md)

## Examples

``` r
# \donttest{
set.seed(1)
X <- matrix(rnorm(500), nrow = 50, ncol = 10)
A <- Matrix::Diagonal(ncol(X))
M <- Matrix::Diagonal(nrow(X))
cv <- blocked_cv_recon_error(X, A, M, k = 3)
# }
```
