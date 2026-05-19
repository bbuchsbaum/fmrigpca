# Gavish–Donoho-style rank heuristic on whitened singular values

Estimates noise level and selects `k` by thresholding the singular
values `s` of the whitened matrix at
`tau = sigma_hat (sqrt(m) + sqrt(n))`.

## Usage

``` r
choose_rank_gd(W)
```

## Arguments

- W:

  Whitened matrix (T x V) or any numeric matrix.

## Value

A list with `k`, `threshold`, `sigma_hat`, and the spectrum `s`.

## See also

[`choose_rank_pa()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_pa.md)
for parallel analysis approach,
[`blocked_cv_recon_error()`](https://bbuchsbaum.github.io/fmrigpca/reference/blocked_cv_recon_error.md)
for cross-validation approach

Other rank selection:
[`blocked_cv_recon_error()`](https://bbuchsbaum.github.io/fmrigpca/reference/blocked_cv_recon_error.md),
[`choose_rank_pa()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_pa.md)

## Examples

``` r
# \donttest{
set.seed(1)
W <- matrix(rnorm(200), nrow = 20, ncol = 10)
gd <- choose_rank_gd(W)
# }
```
