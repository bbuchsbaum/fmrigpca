# Parallel Analysis via time-wise permutation

Computes a null distribution for the leading singular value by permuting
time indices independently within each column and selecting `k` where
observed singular values exceed the `(1 - alpha)` quantile of the null.

## Usage

``` r
choose_rank_pa(W, B = 100L, alpha = 0.05)
```

## Arguments

- W:

  Whitened matrix (T x V).

- B:

  Integer number of permutations (default 100).

- alpha:

  Significance level (default 0.05).

## Value

A list with `k`, `threshold` (null quantile), `null_s1`, and `s`.

## See also

[`choose_rank_gd()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_gd.md)
for Gavish-Donoho approach,
[`blocked_cv_recon_error()`](https://bbuchsbaum.github.io/fmrigpca/reference/blocked_cv_recon_error.md)
for cross-validation approach

Other rank selection:
[`blocked_cv_recon_error()`](https://bbuchsbaum.github.io/fmrigpca/reference/blocked_cv_recon_error.md),
[`choose_rank_gd()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_gd.md)

## Examples

``` r
# \donttest{
set.seed(1)
W <- matrix(rnorm(200), nrow = 20, ncol = 10)
pa <- choose_rank_pa(W, B = 20, alpha = 0.05)
# }
```
