library(testthat)
library(fmrigpca)
library(neuroim2)

test_that("estimate_ar_whitener estimates AR parameters", {
  # Create test data
  test_data <- create_minimal_test_data(dims = c(8, 8, 4), n_time = 50, n_runs = 1)
  nv_run <- test_data$data[[1]]
  
  # Estimate AR(1) whitener
  ar_result <- estimate_ar_whitener(nv_run, wm_mask = test_data$wm, p = 1L)
  
  # Check output structure
  expect_type(ar_result, "list")
  expect_true(all(c("Q", "Sigma", "phi") %in% names(ar_result)))
  
  # Check dimensions
  n_time <- dim(nv_run)[4]
  expect_equal(nrow(ar_result$Q), n_time)
  expect_equal(ncol(ar_result$Q), n_time)
  expect_equal(nrow(ar_result$Sigma), n_time)
  expect_equal(ncol(ar_result$Sigma), n_time)
  
  # Check phi length
  expect_equal(length(ar_result$phi), 1L)  # AR(1)
  
  # phi should be in reasonable range for AR(1)
  expect_true(all(abs(ar_result$phi) < 1))
  
  # Sigma should be positive definite
  expect_no_error(chol(ar_result$Sigma))
})

test_that("estimate_ar_whitener handles higher order AR", {
  test_data <- create_minimal_test_data(dims = c(6, 6, 4), n_time = 60, n_runs = 1)
  nv_run <- test_data$data[[1]]
  
  # AR(2)
  ar2_result <- estimate_ar_whitener(nv_run, p = 2L)
  expect_equal(length(ar2_result$phi), 2L)
  
  # AR(3)
  ar3_result <- estimate_ar_whitener(nv_run, p = 3L)
  expect_equal(length(ar3_result$phi), 3L)
})

test_that("make_temporal_penalty creates valid penalty matrix", {
  n_time <- 30
  
  # Create temporal penalty
  D <- make_temporal_penalty(n_time, lambda_t = 0.3)
  
  # Check dimensions
  expect_equal(nrow(D), n_time)
  expect_equal(ncol(D), n_time)
  
  # Should be symmetric
  expect_true(Matrix::isSymmetric(D))
  
  # Should be positive semi-definite
  eigenvals <- eigen(as.matrix(D), only.values = TRUE)$values
  expect_true(all(eigenvals >= -1e-10))  # Allow small numerical errors
})

test_that("make_temporal_penalty handles different lambda values", {
  n_time <- 20
  
  D1 <- make_temporal_penalty(n_time, lambda_t = 0.1)
  D2 <- make_temporal_penalty(n_time, lambda_t = 0.9)
  
  # Both should be valid
  expect_equal(dim(D1), c(n_time, n_time))
  expect_equal(dim(D2), c(n_time, n_time))
  
  # Higher lambda should mean stronger penalty
  expect_gt(sum(abs(D2)), sum(abs(D1)))
})

test_that("make_frame_weights creates valid weight matrix", {
  n_time <- 40
  
  # Create simple FD and DVARS
  fd <- runif(n_time, 0, 0.5)
  dvars <- runif(n_time, 0.5, 1.5)
  
  # Create frame weights
  W <- make_frame_weights(fd, dvars, fd_thresh = 0.3, dvars_z = 1.0)
  
  # Check dimensions
  expect_equal(nrow(W), n_time)
  expect_equal(ncol(W), n_time)
  
  # Should be diagonal
  expect_true(Matrix::isDiagonal(W))
  
  # Weights should be in (0, 1]
  diag_vals <- Matrix::diag(W)
  expect_true(all(diag_vals > 0))
  expect_true(all(diag_vals <= 1))
})

test_that("make_frame_weights handles missing motion parameters", {
  n_time <- 30
  
  # Only FD
  fd <- runif(n_time, 0, 0.5)
  W_fd <- make_frame_weights(fd, NULL)
  expect_equal(dim(W_fd), c(n_time, n_time))
  
  # Only DVARS
  dvars <- runif(n_time, 0.5, 1.5)
  W_dvars <- make_frame_weights(NULL, dvars)
  expect_equal(dim(W_dvars), c(n_time, n_time))
  
  # Neither (should return identity)
  W_none <- make_frame_weights(NULL, NULL)
  expect_equal(dim(W_none), c(1, 1))  # Default when both NULL
})

test_that("build_temporal_metric creates complete temporal metric", {
  test_data <- create_minimal_test_data(dims = c(8, 8, 4), n_time = 40, n_runs = 1)
  nv_run <- test_data$data[[1]]
  
  # Create motion parameters
  motion <- create_motion_params(40)
  
  # Create M matrix
  M <- build_temporal_metric(
    nv_run, 
    wm_mask = test_data$wm, 
    p = 1L,
    FD = motion$FD,
    DVARS = motion$DVARS,
    lambda_t = 0.3,
    ridge = 1e-6
  )
  
  # Check dimensions
  n_time <- dim(nv_run)[4]
  expect_equal(nrow(M), n_time)
  expect_equal(ncol(M), n_time)
  
  # Should be symmetric
  expect_true(Matrix::isSymmetric(M))
  
  # Should be positive definite (due to ridge)
  expect_no_error(chol(M))
})

test_that("build_temporal_metric handles edge cases", {
  test_data <- create_minimal_test_data(dims = c(6, 6, 4), n_time = 20, n_runs = 1)
  nv_run <- test_data$data[[1]]
  
  # Without WM mask
  M_no_wm <- build_temporal_metric(nv_run, wm_mask = NULL, p = 1L)
  expect_equal(dim(M_no_wm), c(20, 20))
  
  # Without motion parameters
  M_no_motion <- build_temporal_metric(nv_run, p = 1L, FD = NULL, DVARS = NULL)
  expect_equal(dim(M_no_motion), c(20, 20))
  
  # With p = 0 (no AR modeling)
  M_no_ar <- build_temporal_metric(nv_run, p = 0L)
  expect_equal(dim(M_no_ar), c(20, 20))
})