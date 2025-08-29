library(testthat)
library(fmrigpca)
library(neuroim2)

test_that("compute_tsnr calculates temporal SNR correctly", {
  # Create test data
  test_data <- create_minimal_test_data(dims = c(8, 8, 4), n_time = 30, n_runs = 2)
  
  # Compute tSNR
  tsnr <- compute_tsnr(test_data$data)
  
  # Check dimensions
  n_voxels <- nrow(as.matrix(test_data$data[[1]]))
  expect_equal(length(tsnr), n_voxels)
  
  # tSNR should be positive
  expect_true(all(tsnr > 0))
  
  # tSNR should be finite
  expect_true(all(is.finite(tsnr)))
  
  # Reasonable range for tSNR (typically between 1 and 200)
  expect_true(all(tsnr >= 1e-6))  # Lower bound from function
  expect_true(all(tsnr < 1000))   # Reasonable upper bound
})

test_that("compute_tsnr handles single run", {
  test_data <- create_minimal_test_data(dims = c(6, 6, 4), n_time = 20, n_runs = 1)
  
  tsnr <- compute_tsnr(test_data$data)
  
  n_voxels <- nrow(as.matrix(test_data$data[[1]]))
  expect_equal(length(tsnr), n_voxels)
  expect_true(all(tsnr > 0))
})

test_that("make_A_from_vols creates valid column metric", {
  # Create test data
  test_data <- create_minimal_test_data(dims = c(8, 8, 4), n_time = 20)
  
  # Create Laplacian
  L <- make_laplacian(test_data$mask, k = 6, sigma = 2.0, normalized = TRUE)
  
  # Compute tSNR
  tsnr <- compute_tsnr(test_data$data)
  
  # Get mask indices
  mask_idx <- which(as.array(test_data$mask) != 0)
  
  # Create A matrix
  A <- make_A_from_vols(
    test_data$gm, test_data$wm, test_data$csf, L,
    tsnr = tsnr, mask_idx = mask_idx,
    lambda_s = 0.5
  )
  
  # Check dimensions
  n_voxels <- length(mask_idx)
  expect_equal(nrow(A), n_voxels)
  expect_equal(ncol(A), n_voxels)
  
  # A should be symmetric
  expect_true(Matrix::isSymmetric(A))
  
  # A should be positive definite (all eigenvalues > 0)
  # Test with Cholesky decomposition
  expect_no_error(chol(A))
})

test_that("make_A_from_vols handles different tissue weights", {
  test_data <- create_minimal_test_data(dims = c(6, 6, 4))
  L <- make_laplacian(test_data$mask, k = 4)
  mask_idx <- which(as.array(test_data$mask) != 0)
  
  # Test with different alpha, beta, gamma weights
  A1 <- make_A_from_vols(
    test_data$gm, test_data$wm, test_data$csf, L,
    mask_idx = mask_idx,
    alpha = 1.0, beta = 0.5, gamma = 1.0
  )
  
  A2 <- make_A_from_vols(
    test_data$gm, test_data$wm, test_data$csf, L,
    mask_idx = mask_idx,
    alpha = 2.0, beta = 0.1, gamma = 0.5
  )
  
  # Both should be valid
  expect_equal(dim(A1), dim(A2))
  expect_true(Matrix::isSymmetric(A1))
  expect_true(Matrix::isSymmetric(A2))
  
  # But they should be different
  expect_false(all(A1 == A2))
})

test_that("make_A_from_vols handles missing tSNR", {
  test_data <- create_minimal_test_data(dims = c(6, 6, 4))
  L <- make_laplacian(test_data$mask, k = 4)
  mask_idx <- which(as.array(test_data$mask) != 0)
  
  # Without tSNR (should use default of 1)
  A_no_tsnr <- make_A_from_vols(
    test_data$gm, test_data$wm, test_data$csf, L,
    tsnr = NULL, mask_idx = mask_idx
  )
  
  # With tSNR
  tsnr <- compute_tsnr(test_data$data)
  A_with_tsnr <- make_A_from_vols(
    test_data$gm, test_data$wm, test_data$csf, L,
    tsnr = tsnr, mask_idx = mask_idx
  )
  
  # Both should be valid
  expect_equal(dim(A_no_tsnr), dim(A_with_tsnr))
  expect_true(Matrix::isSymmetric(A_no_tsnr))
  expect_true(Matrix::isSymmetric(A_with_tsnr))
})

test_that("make_A_from_vols validates tissue map consistency", {
  test_data <- create_minimal_test_data(dims = c(6, 6, 4))
  
  # Create tissue maps with different space
  wrong_space <- NeuroSpace(c(8, 8, 4), c(2, 2, 2))
  wrong_gm <- NeuroVol(array(0, c(8, 8, 4)), wrong_space)
  
  L <- make_laplacian(test_data$mask, k = 4)
  
  # Should error due to space mismatch
  expect_error(
    make_A_from_vols(wrong_gm, test_data$wm, test_data$csf, L),
    "same NeuroSpace"
  )
})