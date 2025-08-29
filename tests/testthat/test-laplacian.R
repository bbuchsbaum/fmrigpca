library(testthat)
library(fmrigpca)
library(neuroim2)

test_that("make_laplacian creates valid graph Laplacian", {
  # Create small test mask
  mask <- create_spherical_mask(c(10, 10, 6))
  
  # Test basic Laplacian creation
  L <- make_laplacian(mask, k = 6, sigma = 2.0, normalized = TRUE)
  
  # Check dimensions
  n_voxels <- sum(as.array(mask) != 0)
  expect_equal(nrow(L), n_voxels)
  expect_equal(ncol(L), n_voxels)
  
  # Check symmetry
  expect_true(Matrix::isSymmetric(L))
  
  # Check that it's a valid Laplacian (row sums should be ~0 for unnormalized)
  L_unnorm <- make_laplacian(mask, k = 6, sigma = 2.0, normalized = FALSE)
  row_sums <- Matrix::rowSums(L_unnorm)
  expect_true(all(abs(row_sums) < 1e-10))
})

test_that("make_laplacian handles different k values", {
  mask <- create_spherical_mask(c(10, 10, 6))  # Bigger mask for more voxels
  
  # Test with different k values
  L_k3 <- make_laplacian(mask, k = 3, normalized = TRUE)
  L_k6 <- make_laplacian(mask, k = 6, normalized = TRUE)
  
  # More neighbors should lead to denser matrix
  expect_lt(Matrix::nnzero(L_k3), Matrix::nnzero(L_k6))
})

test_that("make_laplacian handles edge cases", {
  # Single voxel mask
  single_vox <- NeuroVol(array(c(0,1,0,0), c(2,2,1)), NeuroSpace(c(2,2,1), c(1,1,1)))
  L_single <- make_laplacian(single_vox, k = 1)
  expect_equal(dim(L_single), c(1, 1))
  
  # Empty mask should error
  empty_mask <- NeuroVol(array(0, c(5,5,5)), NeuroSpace(c(5,5,5), c(1,1,1)))
  expect_error(make_laplacian(empty_mask, k = 3))
})

test_that("make_laplacian requires neighborweights package", {
  skip_if_not_installed("neighborweights")
  
  mask <- create_spherical_mask(c(8, 8, 4))
  n_voxels <- sum(as.array(mask) != 0)
  
  # Should work with neighborweights
  L <- make_laplacian(mask, k = 6, sigma = 2.0, normalized = TRUE)
  expect_equal(nrow(L), n_voxels)
  expect_equal(ncol(L), n_voxels)
  expect_true(Matrix::isSymmetric(L))
})

test_that("normalized vs unnormalized Laplacian properties", {
  mask <- create_spherical_mask(c(10, 10, 6))
  
  L_norm <- make_laplacian(mask, k = 6, normalized = TRUE)
  L_unnorm <- make_laplacian(mask, k = 6, normalized = FALSE)
  
  # Normalized Laplacian eigenvalues should be in [0, 2]
  # (expensive to compute, so only test structure here)
  expect_true(all(Matrix::diag(L_norm) >= 0))
  
  # Unnormalized Laplacian should have zero row sums
  row_sums <- Matrix::rowSums(L_unnorm)
  expect_true(all(abs(row_sums) < 1e-10))
})
