# Combine per-run genpca fits via meta-PCA (MFA)

Aggregates multiple `genpca` fits (e.g., per-run, per-session, or
per-subject) into a single consensus component space using
[`subpca::metapca`](https://rdrr.io/pkg/subpca/man/metapca.html) with
Multiple Factor Analysis (MFA) weighting. This is a “PCA of PCAs” that
operates on the loadings from each individual fit, balancing
contributions across inputs.

## Usage

``` r
fit_subject_metapca(fits, k = 20)
```

## Arguments

- fits:

  List of `genpca` fit objects to combine.

- k:

  Number of meta-components to return (consensus dimensionality).

## Value

A `metapca` object with combined loadings `v` and related metadata.

## Details

Two common workflows are supported in fmrigpca:

1.  Joint fit across runs: stack runs in time, build a block-diagonal
    row metric `M`, and fit once via
    [`fit_subject_genpca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca.md).
    This shares a single spatial column metric `A` and exploits all
    timepoints jointly.

2.  Meta-PCA across individual fits: fit each run (or session/subject)
    separately to obtain `genpca` objects, then combine them with
    meta-PCA.

Meta-PCA uses MFA weighting to avoid any single fit dominating due to
scale or sample size differences, and returns a `metapca` object with
combined loadings (accessible as `$v`). Use this when:

- Runs/sessions are heterogeneous and you prefer to fit them
  independently (e.g., different preprocessing, durations, or
  acquisition settings)

- Combining across subjects after per-subject genpca fits

Assumptions and tips:

- Loadings should be defined on the same voxel/parcel index set (e.g.,
  same mask and order). If they differ, align to a common subset before
  fitting.

- Keep `k` small and comparable across inputs; `metapca` selects a
  consensus space up to `ncomp = k`.

- See the subpca package for algorithmic details of MFA combination.

## See also

[`fit_subject_genpca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca.md)
and
[`fit_subject_genpca_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca_parcel.md)
for generating individual fits to combine

Other main fitting:
[`fit_subject_genpca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca.md),
[`fit_subject_genpca_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca_parcel.md)

## Examples

``` r
# \donttest{
# Per-run fits (each list contains a single run)
fits <- lapply(seq_along(nv_list), function(r) {
  fit_subject_genpca(
    nv_list = list(nv_list[[r]]),
    mask_vol = mask,
    gm_vol = gm, wm_vol = wm, csf_vol = csf,
    k = 5, p_ar = 1L
  )
})
#> Error: object 'nv_list' not found
# Combine via meta-PCA (MFA weighting)
meta <- fit_subject_metapca(fits, k = 3)
#> Error: object 'fits' not found
str(meta$v)  # consensus loadings (voxels x k)
#> Error: object 'meta' not found
# }
```
