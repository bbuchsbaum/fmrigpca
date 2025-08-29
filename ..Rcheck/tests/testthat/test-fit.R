library(testthat)
library(fmrigpca)
library(neuroim2)
library(genpca)

test_that("fit_subject_genpca runs complete pipeline", {
  skip_if_not_installed("genpca")
  
  # Create small test dataset with enough voxels
  test_data <- create_minimal_test_data(dims = c(10, 10, 6), n_time = 30, n_runs = 2)
  
  # Run genpca fitting
  result <- fit_subject_genpca(
    nv_list = test_data$data,
    mask_vol = test_data$mask,
    gm_vol = test_data$gm,
    wm_vol = test_data$wm,
    csf_vol = test_data$csf,
    FD_list = test_data$FD,
    DVARS_list = test_data$DVARS,
    k = 5,  # Small number of components for testing
    p_ar = 1L,
    lambda_s = 0.5,
    lambda_t = 0.3,
    lap_k = 3,
    lap_sigma = 2.0
  )
  
  # Check that result is a genpca object
  expect_s3_class(result, "genpca")
  
  # Check dimensions
  expect_equal(ncol(result$v), 5)  # k components
  n_voxels <- sum(as.array(test_data$mask) != 0)
  expect_equal(nrow(result$v), n_voxels)
  
  # Check that scores exist
  expect_true(!is.null(result$scores))
  total_time <- sum(sapply(test_data$data, function(x) dim(x)[4]))
  expect_equal(nrow(result$scores), total_time)
  expect_equal(ncol(result$scores), 5)
})

test_that("fit_subject_genpca validates inputs", {
  test_data <- create_minimal_test_data(dims = c(8, 8, 5), n_time = 20, n_runs = 2)
  
  # Different space for tissue maps should error
  wrong_space <- NeuroSpace(c(8, 8, 4), c(2, 2, 2))
  wrong_gm <- NeuroVol(array(0, c(8, 8, 4)), wrong_space)
  
  expect_error(
    fit_subject_genpca(
      nv_list = test_data$data,
      mask_vol = test_data$mask,
      gm_vol = wrong_gm,
      wm_vol = test_data$wm,
      csf_vol = test_data$csf,
      k = 3
    ),
    "same NeuroSpace"
  )
  
  # Empty mask should error
  empty_mask <- NeuroVol(array(0, dim(test_data$mask)), space(test_data$mask))
  expect_error(
    fit_subject_genpca(
      nv_list = test_data$data,
      mask_vol = empty_mask,
      gm_vol = test_data$gm,
      wm_vol = test_data$wm,
      csf_vol = test_data$csf,
      k = 3
    ),
    "Mask has no voxels"
  )
})

test_that("fit_subject_genpca handles single run", {
  skip_if_not_installed("genpca")
  
  test_data <- create_minimal_test_data(dims = c(8, 8, 5), n_time = 25, n_runs = 1)
  
  result <- fit_subject_genpca(
    nv_list = test_data$data,
    mask_vol = test_data$mask,
    gm_vol = test_data$gm,
    wm_vol = test_data$wm,
    csf_vol = test_data$csf,
    k = 3,
    p_ar = 1L
  )
  
  expect_s3_class(result, "genpca")
  expect_equal(ncol(result$v), 3)
})

test_that("fit_subject_genpca works without motion parameters", {
  skip_if_not_installed("genpca")
  
  test_data <- create_minimal_test_data(dims = c(8, 8, 5), n_time = 20, n_runs = 2)
  
  # Run without FD and DVARS
  result <- fit_subject_genpca(
    nv_list = test_data$data,
    mask_vol = test_data$mask,
    gm_vol = test_data$gm,
    wm_vol = test_data$wm,
    csf_vol = test_data$csf,
    FD_list = NULL,
    DVARS_list = NULL,
    k = 3
  )
  
  expect_s3_class(result, "genpca")
  expect_equal(ncol(result$v), 3)
})

test_that("fit_subject_metapca combines multiple fits", {
  skip_if_not_installed("genpca")
  skip_if_not_installed("subpca")
  
  test_data <- create_minimal_test_data(dims = c(8, 8, 5), n_time = 20, n_runs = 3)
  
  # Create individual fits for each run
  fits <- list()
  for (r in 1:3) {
    fits[[r]] <- fit_subject_genpca(
      nv_list = list(test_data$data[[r]]),  # Single run
      mask_vol = test_data$mask,
      gm_vol = test_data$gm,
      wm_vol = test_data$wm,
      csf_vol = test_data$csf,
      k = 4
    )
  }
  
  # Combine with metapca
  meta_result <- fit_subject_metapca(fits, k = 3)
  
  # Check result
  expect_s3_class(meta_result, "metapca")
  expect_equal(ncol(meta_result$v), 3)
})

test_that("fit_subject_genpca handles different parameters", {
  skip_if_not_installed("genpca")
  
  test_data <- create_minimal_test_data(dims = c(8, 8, 5), n_time = 20, n_runs = 1)
  
  # Test with different lambda values
  result1 <- fit_subject_genpca(
    nv_list = test_data$data,
    mask_vol = test_data$mask,
    gm_vol = test_data$gm,
    wm_vol = test_data$wm,
    csf_vol = test_data$csf,
    k = 3,
    lambda_s = 0.1,
    lambda_t = 0.1
  )
  
  result2 <- fit_subject_genpca(
    nv_list = test_data$data,
    mask_vol = test_data$mask,
    gm_vol = test_data$gm,
    wm_vol = test_data$wm,
    csf_vol = test_data$csf,
    k = 3,
    lambda_s = 0.9,
    lambda_t = 0.9
  )
  
  # Both should work
  expect_s3_class(result1, "genpca")
  expect_s3_class(result2, "genpca")
  
  # Results should be different due to different regularization
  expect_false(all(result1$v == result2$v))
})