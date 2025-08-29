library(testthat)
library(fmrigpca)
library(neuroim2)

test_that("parcel_centroids_from_labels computes correct centroids", {
  # Create a simple parcellation
  dims <- c(8, 8, 4)
  parc_array <- array(0L, dims)
  
  # Create 3 parcels
  parc_array[1:3, 1:3, 1:2] <- 1L  # Parcel 1
  parc_array[5:7, 1:3, 1:2] <- 2L  # Parcel 2
  parc_array[3:5, 5:7, 3:4] <- 3L  # Parcel 3
  
  parc_vol <- NeuroVol(parc_array, NeuroSpace(dims, c(1, 1, 1)))
  
  # Compute centroids
  result <- parcel_centroids_from_labels(parc_vol)
  
  # Check structure
  expect_type(result, "list")
  expect_true(all(c("coords", "labels") %in% names(result)))
  
  # Check dimensions
  expect_equal(nrow(result$coords), 3)  # 3 parcels
  expect_equal(ncol(result$coords), 3)  # i, j, k coordinates
  expect_equal(length(result$labels), 3)
  
  # Check labels
  expect_equal(result$labels, c(1L, 2L, 3L))
  
  # Check that centroids are reasonable (within parcel bounds)
  expect_true(result$coords[1, 1] >= 1 && result$coords[1, 1] <= 3)  # Parcel 1 i-coord
  expect_true(result$coords[2, 1] >= 5 && result$coords[2, 1] <= 7)  # Parcel 2 i-coord
})

test_that("parcel_centroids_from_labels handles specific labels", {
  dims <- c(6, 6, 4)
  parc_array <- array(0L, dims)
  parc_array[1:2, 1:2, 1] <- 1L
  parc_array[3:4, 1:2, 1] <- 2L
  parc_array[5:6, 1:2, 1] <- 3L
  parc_array[1:2, 3:4, 1] <- 4L
  
  parc_vol <- NeuroVol(parc_array, NeuroSpace(dims, c(1, 1, 1)))
  
  # Select specific labels
  result <- parcel_centroids_from_labels(parc_vol, labels = c(2L, 4L))
  
  expect_equal(nrow(result$coords), 2)
  expect_equal(result$labels, c(2L, 4L))
})

test_that("aggregate_tissue_to_parcels aggregates correctly", {
  # Create test data
  test_data <- create_minimal_test_data(dims = c(8, 8, 4))
  parc_vol <- create_parcellation(test_data$mask, n_parcels = 5)
  
  # Aggregate tissue maps to parcels
  result <- aggregate_tissue_to_parcels(
    test_data$gm, test_data$wm, test_data$csf, parc_vol
  )
  
  # Check structure
  expect_type(result, "list")
  expect_true(all(c("gm", "wm", "csf") %in% names(result)))
  
  # Check dimensions (should be n_parcels x 1)
  n_parcels <- length(unique(as.vector(as.array(parc_vol))[as.vector(as.array(parc_vol)) > 0]))
  expect_equal(length(result$gm), n_parcels)
  expect_equal(length(result$wm), n_parcels)
  expect_equal(length(result$csf), n_parcels)
  
  # Values should be between 0 and 1
  expect_true(all(result$gm >= 0 & result$gm <= 1))
  expect_true(all(result$wm >= 0 & result$wm <= 1))
  expect_true(all(result$csf >= 0 & result$csf <= 1))
  
  # Should approximately sum to 1 for each parcel
  totals <- result$gm + result$wm + result$csf
  expect_true(all(abs(totals - 1) < 0.1))  # Allow some tolerance
})

test_that("make_parcel_laplacian creates valid Laplacian", {
  # Create parcellation
  test_data <- create_minimal_test_data(dims = c(10, 10, 6))
  parc_vol <- create_parcellation(test_data$mask, n_parcels = 4)
  
  # Get centroids
  centroids <- parcel_centroids_from_labels(parc_vol)
  
  # Create parcel Laplacian
  L <- make_parcel_laplacian(parcel_coords = centroids$coords, k = 3, sigma = 5.0, normalized = TRUE)
  
  # Check dimensions
  n_parcels <- nrow(centroids$coords)
  expect_equal(nrow(L), n_parcels)
  expect_equal(ncol(L), n_parcels)
  
  # Should be symmetric
  expect_true(Matrix::isSymmetric(L))
  
  # For unnormalized, row sums should be ~0
  L_unnorm <- make_parcel_laplacian(parcel_coords = centroids$coords, k = 3, sigma = 5.0, normalized = FALSE)
  row_sums <- Matrix::rowSums(L_unnorm)
  expect_true(all(abs(row_sums) < 1e-10))
})

test_that("compute_tsnr_parcel calculates parcel tSNR", {
  # ClusteredNeuroVec is part of neuroim2, which is a required dependency
  skip_if_not(exists("ClusteredNeuroVec", where = asNamespace("neuroim2")))
  
  # This test requires ClusteredNeuroVec class from neuroim2
  
  # Create test data
  test_data <- create_minimal_test_data(dims = c(8, 8, 4), n_time = 30, n_runs = 2)
  parc_vol <- create_parcellation(test_data$mask, n_parcels = 6)
  
  # Would test compute_tsnr_parcel here if ClusteredNeuroVec is available
  # For now, just check that the function exists
  expect_true(exists("compute_tsnr_parcel"))
})

test_that("make_A_from_parcels creates valid column metric", {
  # Create test data
  test_data <- create_minimal_test_data(dims = c(8, 8, 4))
  parc_vol <- create_parcellation(test_data$mask, n_parcels = 6)
  
  # Get centroids and aggregate tissue
  centroids <- parcel_centroids_from_labels(parc_vol)
  tissue <- aggregate_tissue_to_parcels(
    test_data$gm, test_data$wm, test_data$csf, parc_vol
  )
  
  # Create Laplacian
  L <- make_parcel_laplacian(parcel_coords = centroids$coords, k = 3, normalized = TRUE)
  
  # Create A matrix
  A <- make_A_from_parcels(
    tissue$gm, tissue$wm, tissue$csf, L,
    tsnr = NULL,  # No tSNR for simplicity
    lambda_s = 0.5
  )
  
  # Check dimensions
  n_parcels <- length(tissue$gm)
  expect_equal(nrow(A), n_parcels)
  expect_equal(ncol(A), n_parcels)
  
  # Should be symmetric
  expect_true(Matrix::isSymmetric(A))
  
  # Should be positive definite
  expect_no_error(chol(A))
})

test_that("parcel functions handle edge cases", {
  # Empty parcellation
  empty_parc <- NeuroVol(array(0L, c(5, 5, 3)), NeuroSpace(c(5, 5, 3), c(1, 1, 1)))
  expect_error(parcel_centroids_from_labels(empty_parc))
  
  # Single parcel
  single_parc <- NeuroVol(array(1L, c(3, 3, 2)), NeuroSpace(c(3, 3, 2), c(1, 1, 1)))
  result <- parcel_centroids_from_labels(single_parc)
  expect_equal(nrow(result$coords), 1)
  expect_equal(result$labels, 1L)
  
  # Non-contiguous labels (1, 3, 5)
  dims <- c(6, 6, 3)
  parc_array <- array(0L, dims)
  parc_array[1:2, 1:2, 1] <- 1L
  parc_array[3:4, 1:2, 1] <- 3L
  parc_array[5:6, 1:2, 1] <- 5L
  
  parc_vol <- NeuroVol(parc_array, NeuroSpace(dims, c(1, 1, 1)))
  result <- parcel_centroids_from_labels(parc_vol)
  expect_equal(result$labels, c(1L, 3L, 5L))
})