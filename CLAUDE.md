# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## Project Overview

`fmrigpca` is an R package for design-free fMRI dimensionality reduction
using generalized PCA (genpca). It provides tools to build row/column
metrics for genpca on fMRI data using NeuroVec/NeuroVol from the
neuroim2 package.

## Key Commands

### Building and Installing

``` r

# Install package dependencies from GitHub
devtools::install_github("bbuchsbaum/genpca")
devtools::install_github("bbuchsbaum/neuroim2")
devtools::install_github("bbuchsbaum/graphweights")
devtools::install_github("bbuchsbaum/subpca")

# Install package locally with vignettes
devtools::install_local("fmrigpca_0.1.3.tar.gz", build_vignettes = TRUE, upgrade = "never")

# Build documentation
devtools::document()

# Check package
devtools::check()

# Build vignettes
devtools::build_vignettes()
```

### Running Tests

``` r

# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-fit.R")
testthat::test_file("tests/testthat/test-laplacian.R")
testthat::test_file("tests/testthat/test-parcel.R")
testthat::test_file("tests/testthat/test-rank.R")
testthat::test_file("tests/testthat/test-spatial.R")
testthat::test_file("tests/testthat/test-temporal.R")
```

## Architecture

The package is organized around two main data processing paradigms:

### Voxel-wise Processing (NeuroVec)

- **Core functions**:
  [`fit_subject_genpca()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca.md) -
  main entry point for voxel-wise analysis
- **Spatial metrics**:
  [`make_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_laplacian.md)
  creates graph Laplacian from mask,
  [`build_spatial_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_spatial_metric.md)
  combines tissue maps with Laplacian
- **Temporal metrics**:
  [`estimate_ar_whitener()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener.md)
  for AR(p) whitening,
  [`build_temporal_metric()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric.md)
  creates temporal penalty matrices
- **Key concept**: Works directly with voxel-level data in `NeuroVec`
  format

### Parcel-wise Processing (ClusteredNeuroVec)

- **Core functions**:
  [`fit_subject_genpca_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/fit_subject_genpca_parcel.md) -
  main entry point for parcel-wise analysis
- **Parcel utilities**:
  [`parcel_centroids_from_labels()`](https://bbuchsbaum.github.io/fmrigpca/reference/parcel_centroids_from_labels.md),
  [`aggregate_tissue_to_parcels()`](https://bbuchsbaum.github.io/fmrigpca/reference/aggregate_tissue_to_parcels.md)
- **Spatial metrics**:
  [`make_parcel_laplacian()`](https://bbuchsbaum.github.io/fmrigpca/reference/make_parcel_laplacian.md)
  creates parcel-level Laplacian
- **Temporal metrics**:
  [`estimate_ar_whitener_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/estimate_ar_whitener_parcel.md),
  [`build_temporal_metric_parcel()`](https://bbuchsbaum.github.io/fmrigpca/reference/build_temporal_metric_parcel.md)
  for parcel-specific whitening
- **Key concept**: Aggregates data to parcels for computational
  efficiency

### Cross-cutting Concerns

- **Rank selection**:
  [`choose_rank_gd()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_gd.md)
  (Gabriel-style cross-validation),
  [`choose_rank_pa()`](https://bbuchsbaum.github.io/fmrigpca/reference/choose_rank_pa.md)
  (parallel analysis)
- **Stability analysis**:
  [`subspace_principal_angles()`](https://bbuchsbaum.github.io/fmrigpca/reference/subspace_principal_angles.md),
  [`procrustes_distance()`](https://bbuchsbaum.github.io/fmrigpca/reference/procrustes_distance.md)
  for comparing subspaces
- **Integration**: Both voxel and parcel approaches feed into
  [`genpca::genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.html)
  with custom A (spatial) and M (temporal) matrices

## Key Dependencies

The package relies heavily on: - **neuroim2**: Provides
`NeuroVec`/`NeuroVol` neuroimaging data structures - **genpca**: Core
generalized PCA implementation - **Matrix**: Sparse matrix operations
for efficiency - **neighborweights**: Graph weight calculations -
**subpca** (optional): For meta-PCA combination of multiple fits

## Development Notes

- The package uses S3 methods and exports functions via NAMESPACE (no S4
  or R6 classes)
- Test suite exists in `tests/testthat/` with comprehensive coverage
  across all modules
- Vignettes in `vignettes/` demonstrate typical workflows
- Test files cover: fit.R, laplacian.R, parcel.R, rank.R, spatial.R, and
  temporal.R modules
