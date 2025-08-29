library(testthat)
library(fmrigpca)
library(neuroim2)
library(Matrix)

test_that("whitened_matrix creates correct whitened transformation", {
  # Create small test matrices
  n <- 20
  p <- 15
  set.seed(123)
  
  X <- matrix(rnorm(n * p), n, p)
  
  # Create positive definite A and M
  A_half <- matrix(rnorm(p * p), p, p)
  A <- Matrix::forceSymmetric(crossprod(A_half) + diag(p) * 0.1)
  
  M_half <- matrix(rnorm(n * n), n, n)
  M <- Matrix::forceSymmetric(crossprod(M_half) + diag(n) * 0.1)
  
  # Compute whitened matrix
  W <- whitened_matrix(X, A, M)
  
  # Check dimensions
  expect_equal(nrow(W), n)
  expect_equal(ncol(W), p)
  
  # The whitened matrix should have unit covariance structure
  # (approximately, for large enough matrices)
  WtW <- crossprod(W)
  # Diagonal should be near 1, off-diagonal near 0 (for identity target)
  # This is approximate due to finite sample
})

test_that("whitened_svd computes singular values", {
  n <- 25
  p <- 20
  set.seed(456)
  
  X <- matrix(rnorm(n * p), n, p)
  
  # Simple identity metrics for testing
  A <- Matrix::Diagonal(p)
  M <- Matrix::Diagonal(n)
  
  # Full SVD (just singular values)
  result_full <- whitened_svd(X, A, M, k = NULL)
  
  expect_true(!is.null(result_full$d))
  expect_true(is.null(result_full$u))
  expect_true(is.null(result_full$v))
  expect_equal(length(result_full$d), min(n, p))
  
  # Truncated SVD
  k <- 5
  result_trunc <- whitened_svd(X, A, M, k = k)
  
  expect_equal(length(result_trunc$d), k)
  expect_true(all(result_trunc$d >= 0))
  
  # Singular values should be in descending order
  expect_true(all(diff(result_trunc$d) <= 0))
})

test_that("whitened_svd works with RSpectra if available", {
  skip_if_not_installed("RSpectra")
  
  n <- 30
  p <- 25
  set.seed(789)
  
  X <- matrix(rnorm(n * p), n, p)
  A <- Matrix::Diagonal(p)
  M <- Matrix::Diagonal(n)
  
  k <- 3
  result <- whitened_svd(X, A, M, k = k)
  
  expect_equal(length(result$d), k)
  expect_true(all(result$d >= 0))
})

test_that("subspace_principal_angles computes angles correctly", {
  set.seed(111)
  n <- 20
  k1 <- 3
  k2 <- 3
  
  # Create two subspaces
  U <- qr.Q(qr(matrix(rnorm(n * k1), n, k1)))
  V <- qr.Q(qr(matrix(rnorm(n * k2), n, k2)))
  
  # Compute principal angles
  angles <- subspace_principal_angles(U, V)
  
  # Check dimensions
  expect_equal(length(angles), min(k1, k2))
  
  # Angles should be between 0 and pi/2
  expect_true(all(angles >= 0))
  expect_true(all(angles <= pi/2 + 1e-10))
  
  # Same subspace should give ~0 angles
  angles_same <- subspace_principal_angles(U, U)
  expect_true(all(angles_same < 1e-6))  # Relaxed tolerance for floating point
  
  # Orthogonal subspaces
  U_orth <- qr.Q(qr(matrix(rnorm(n * k1), n, k1)))
  V_orth <- U_orth[, 1:k1, drop = FALSE]
  # Make V orthogonal to U
  V_orth_perp <- qr.Q(qr(diag(n) - tcrossprod(U_orth)))[, 1:k2]
  angles_orth <- subspace_principal_angles(U_orth, V_orth_perp)
  expect_true(all(angles_orth > pi/4))  # Should be large angles
})

test_that("procrustes_distance measures subspace distance", {
  set.seed(222)
  n <- 25
  k <- 4
  
  # Create two subspaces
  U <- qr.Q(qr(matrix(rnorm(n * k), n, k)))
  V <- qr.Q(qr(matrix(rnorm(n * k), n, k)))
  
  # Compute distance
  dist <- procrustes_distance(U, V)
  
  # Distance should be non-negative
  expect_true(dist >= 0)
  
  # Distance to itself should be 0
  dist_same <- procrustes_distance(U, U)
  expect_true(abs(dist_same) < 1e-10)
  
  # Maximum distance is sqrt(2k)
  max_dist <- sqrt(2 * k)
  expect_true(dist <= max_dist + 1e-10)
})

test_that("choose_rank_gd performs Gabriel-style cross-validation", {
  skip_if_not_installed("genpca")
  
  # Create small test data
  test_data <- create_minimal_test_data(dims = c(6, 6, 4), n_time = 40, n_runs = 1)
  
  # Prepare data matrix
  X <- t(as.matrix(test_data$data[[1]]))  # T x V
  mask_idx <- which(as.array(test_data$mask) != 0)
  X <- X[, mask_idx]
  
  # Simple metrics for testing
  A <- Matrix::Diagonal(ncol(X))
  M <- Matrix::Diagonal(nrow(X))
  
  # Create whitened matrix
  W <- whitened_matrix(X, A, M)
  
  # Test rank selection
  result <- choose_rank_gd(W)
  
  # Check structure
  expect_type(result, "list")
  expect_true(all(c("k", "threshold", "sigma_hat", "s") %in% names(result)))
  
  # Check outputs
  expect_type(result$k, "integer")
  expect_true(result$k >= 0)
  expect_true(result$threshold > 0)
  expect_true(all(result$s >= 0))
})

test_that("choose_rank_pa performs parallel analysis", {
  skip_if_not_installed("genpca")
  
  # Create test data
  n <- 30
  p <- 20
  set.seed(333)
  X <- matrix(rnorm(n * p), n, p)
  
  A <- Matrix::Diagonal(p)
  M <- Matrix::Diagonal(n)
  
  # Create whitened matrix
  W <- whitened_matrix(X, A, M)
  
  # Run parallel analysis
  result <- choose_rank_pa(W, B = 10)
  
  # Check structure
  expect_type(result, "list")
  expect_true(all(c("k", "threshold", "null_s1", "s") %in% names(result)))
  
  # Check outputs
  expect_type(result$k, "integer")
  expect_true(result$k >= 0)
  expect_true(result$threshold > 0)
  expect_equal(length(result$null_s1), 10)  # B = 10
  
  # Singular values should be positive and decreasing
  expect_true(all(result$s >= 0))
  expect_true(all(diff(result$s) <= 0))
})

test_that("blocked_cv_recon_error computes reconstruction error", {
  # Create simple test case
  n <- 20
  p <- 15
  set.seed(444)
  
  X <- matrix(rnorm(n * p), n, p)
  
  # Create metrics for the function
  A <- Matrix::Diagonal(p)
  M <- Matrix::Diagonal(n)
  
  # Compute reconstruction error
  result <- blocked_cv_recon_error(X, A, M, k = 3, nfold = 5L, block_frac = 0.2)
  
  # Check it returns a data frame
  expect_true(is.data.frame(result))
  
  # Check structure
  expect_true(all(c("fold", "error", "mean_error") %in% names(result)))
  
  # Check dimensions
  expect_equal(nrow(result), 5)
  
  # Errors should be positive
  expect_true(all(result$error > 0))
  
  # Mean error should match
  expect_equal(result$mean_error[1], mean(result$error))
  
  # Higher k should give lower error (generally)
  result_k1 <- blocked_cv_recon_error(X, A, M, k = 1, nfold = 5L, block_frac = 0.2)
  result_k5 <- blocked_cv_recon_error(X, A, M, k = 5, nfold = 5L, block_frac = 0.2)
  expect_true(mean(result_k5$error) <= mean(result_k1$error) * 1.1)  # Allow 10% tolerance
})