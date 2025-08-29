# Helper functions for generating test data

#' Create a spherical brain mask
#' @param dims Dimensions of the volume (default c(16, 16, 8))
#' @param radius_fraction Fraction of smallest dimension for radius (default 0.6)
create_spherical_mask <- function(dims = c(16, 16, 8), radius_fraction = 0.6) {
  mask_array <- array(FALSE, dims)
  center <- (dims + 1) / 2  # Center of the volume
  radius <- min(dims) * radius_fraction / 2
  
  for (i in 1:dims[1]) {
    for (j in 1:dims[2]) {
      for (k in 1:dims[3]) {
        # Simple Euclidean distance
        dist <- sqrt(sum((c(i, j, k) - center)^2))
        if (dist <= radius) {
          mask_array[i,j,k] <- TRUE
        }
      }
    }
  }
  
  neuroim2::NeuroVol(mask_array, neuroim2::NeuroSpace(dims, c(2, 2, 2)))
}

#' Create synthetic tissue probability maps
#' @param mask_vol Brain mask volume
#' @param gm_center_fraction Where GM is centered (0.5 = middle)
create_tissue_maps <- function(mask_vol, gm_center_fraction = 0.5) {
  dims <- dim(mask_vol)
  mask_array <- as.array(mask_vol) != 0
  mask_idx <- which(mask_array)
  
  # Get voxel coordinates
  lin <- mask_idx - 1
  ii <- (lin %% dims[1]) + 1
  jj <- ((lin %/% dims[1]) %% dims[2]) + 1
  kk <- (lin %/% (dims[1] * dims[2])) + 1
  
  # Compute distance from center for each voxel
  center <- dims / 2
  dists <- sqrt(rowSums(sweep(cbind(ii, jj, kk), 2, center, "-")^2))
  max_dist <- max(dists)
  
  # Create tissue maps based on distance
  # GM peaks at intermediate distance
  # WM is central
  # CSF is peripheral
  
  gm_array <- wm_array <- csf_array <- array(0, dims)
  
  if (length(mask_idx) > 0) {
    for (idx in seq_along(mask_idx)) {
      vox <- mask_idx[idx]
      d_norm <- dists[idx] / max_dist
      
      # WM peaks at center
      wm_array[vox] <- exp(-4 * d_norm^2)
      
      # GM peaks at intermediate distance
      gm_array[vox] <- exp(-8 * (d_norm - gm_center_fraction)^2)
      
      # CSF increases with distance
      csf_array[vox] <- pmin(1, d_norm^2)
    }
    
    # Normalize so they sum to ~1 at each voxel
    total <- gm_array + wm_array + csf_array
    total[total == 0] <- 1
    
    gm_array <- gm_array / total
    wm_array <- wm_array / total
    csf_array <- csf_array / total
  }
  
  list(
    gm = neuroim2::NeuroVol(gm_array, neuroim2::space(mask_vol)),
    wm = neuroim2::NeuroVol(wm_array, neuroim2::space(mask_vol)),
    csf = neuroim2::NeuroVol(csf_array, neuroim2::space(mask_vol))
  )
}

#' Create a simple parcellation volume
#' @param mask_vol Brain mask volume
#' @param n_parcels Number of parcels to create
create_parcellation <- function(mask_vol, n_parcels = 10) {
  dims <- dim(mask_vol)
  mask_array <- as.array(mask_vol) != 0
  mask_idx <- which(mask_array)
  n_vox <- length(mask_idx)
  
  if (n_parcels > n_vox) {
    n_parcels <- max(1, floor(n_vox / 2))  # Ensure parcels < voxels for kmeans
  }
  
  # Simple k-means based parcellation
  # Get coordinates
  lin <- mask_idx - 1
  ii <- (lin %% dims[1]) + 1
  jj <- ((lin %/% dims[1]) %% dims[2]) + 1
  kk <- (lin %/% (dims[1] * dims[2])) + 1
  coords <- cbind(ii, jj, kk)
  
  # Perform k-means clustering
  set.seed(123)
  km <- kmeans(coords, centers = n_parcels, nstart = 5)
  
  # Create parcellation volume
  parc_array <- array(0L, dims)
  parc_array[mask_idx] <- km$cluster
  
  neuroim2::NeuroVol(parc_array, neuroim2::space(mask_vol))
}

#' Generate simulated fMRI data with known properties
#' @param mask_vol Brain mask volume
#' @param n_time Number of time points
#' @param n_runs Number of runs to generate
#' @param seed Random seed for reproducibility
create_test_fmri_data <- function(mask_vol, n_time = 50, n_runs = 2, seed = 42) {
  nv_list <- list()
  
  for (r in 1:n_runs) {
    set.seed(seed + r - 1)
    nv_list[[r]] <- neuroim2::simulate_fmri(
      mask = mask_vol,
      n_time = n_time,
      spatial_fwhm = 4,      # Moderate smoothing
      ar_mean = 0.4,         # Typical AR(1)
      ar_sd = 0.08,          # Keep AR coefficients reasonable
      noise_sd = 1.0,
      hetero_strength = 0.2, # Reduced heteroscedasticity to avoid numerical issues
      global_amp = 0.1,      # Small global signal
      n_factors = 3,         # Few latent components
      seed = seed + r - 1,
      return_centered = TRUE
    )
  }
  
  nv_list
}

#' Generate simulated motion parameters
#' @param n_time Number of time points
#' @param mean_fd Mean framewise displacement
#' @param mean_dvars Mean DVARS
create_motion_params <- function(n_time, mean_fd = 0.2, mean_dvars = 1.0) {
  # FD: Framewise displacement (motion)
  fd <- abs(rnorm(n_time, mean = mean_fd, sd = mean_fd/2))
  fd[fd < 0] <- 0
  
  # DVARS: RMS of derivative (signal change)
  dvars <- abs(rnorm(n_time, mean = mean_dvars, sd = mean_dvars/3))
  dvars[dvars < 0] <- 0.1
  
  list(FD = fd, DVARS = dvars)
}

#' Create a minimal test dataset for quick testing
#' @param dims Dimensions for the brain volume
#' @param n_time Number of time points
#' @param n_runs Number of runs
create_minimal_test_data <- function(dims = c(10, 10, 6), n_time = 20, n_runs = 1) {
  # Create mask
  mask_vol <- create_spherical_mask(dims, radius_fraction = 0.7)
  
  # Create tissue maps
  tissue_maps <- create_tissue_maps(mask_vol, gm_center_fraction = 0.5)
  
  # Create simulated fMRI data
  nv_list <- create_test_fmri_data(mask_vol, n_time = n_time, n_runs = n_runs)
  
  # Create motion parameters
  motion_list <- list()
  for (r in 1:n_runs) {
    motion_list[[r]] <- create_motion_params(n_time)
  }
  
  list(
    mask = mask_vol,
    gm = tissue_maps$gm,
    wm = tissue_maps$wm,
    csf = tissue_maps$csf,
    data = nv_list,
    FD = lapply(motion_list, function(x) x$FD),
    DVARS = lapply(motion_list, function(x) x$DVARS)
  )
}