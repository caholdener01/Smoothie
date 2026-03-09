# Copyright (c) 2026 Chase Holdener
# Licensed under the MIT License. See LICENSE file for details.

import pandas as pd
import numpy as np
import anndata as ad
import time
import random
import gc
import sys
import multiprocessing as mp

from scipy.sparse import issparse, csr_matrix, csc_matrix, coo_matrix, lil_matrix
from scipy.spatial import cKDTree
from scipy.stats import pearsonr


## FUNCTION FOR SHUFFLING ADATA IN BINS

# Shuffles adata in bins, keeping the relative position of coords within each bin preserved
def bin_shuffle_adata(adata, bin_width, seed):

    min_x = np.min(adata.obsm['spatial'][:, 0])
    min_y = np.min(adata.obsm['spatial'][:, 1])
    adata.obsm['spatial'][:, 0] -= min_x
    adata.obsm['spatial'][:, 1] -= min_y
    
    # Place each coord into a bin
    adata.obsm['bin_coords'] = np.stack([
        np.ceil(adata.obsm['spatial'][:, 0] / bin_width),
        np.ceil(adata.obsm['spatial'][:, 1] / bin_width)
    ], axis=1)
    
    np.random.seed(seed)
    # Make dictionary for shuffling
    unique_bins, inverse_indices = np.unique(adata.obsm['bin_coords'], axis=0, return_inverse=True)
    shuffled_unique_bins = np.random.permutation(unique_bins)
    
    adata.obsm['shuffled_bin_coords'] = shuffled_unique_bins[inverse_indices]
    
    # Update spatial coords of adata
    adata.obsm['spatial'] = \
        adata.obsm['spatial'] + (bin_width * (adata.obsm['shuffled_bin_coords'] - adata.obsm['bin_coords']))


## FUNCTIONS FOR CONSTRUCTING GRID

# Makes hexagonal grid within range of min and max x and y values
def make_grid(max_X, min_X, max_Y, min_Y, stride):

    # get sqrt(3)/2 * stride
    stride_y_dir = (np.sqrt(3) / 2.0) * np.float32(stride)
    
    x_vals = np.arange(min_X, max_X + stride, stride)
    y_vals = np.arange(min_Y, max_Y + stride_y_dir, stride_y_dir)

    grid_points = []
    # Shift every other row horizontally by stride * 0.5
    for i, row in enumerate(y_vals):
        if i % 2 != 0:
            # shifted row
            grid_points.extend([(x + (stride * 0.5), y) for x, y in zip(x_vals[:-1], [row] * (len(x_vals) - 1))])
        else:
            # non-shifted row
            grid_points.extend([(x, y) for x, y in zip(x_vals[:-1], [row] * (len(x_vals) - 1))])

    return grid_points


# Makes square grid within range of min and max x and y values
def make_square_grid(max_X, min_X, max_Y, min_Y, stride):
    x_vals = np.arange(min_X, max_X + stride, stride)
    y_vals = np.arange(min_Y, max_Y + stride, stride)

    X, Y = np.meshgrid(x_vals[:-1], y_vals[:-1])  # grid of lower-left corners
    grid_points = np.column_stack((X.ravel(), Y.ravel()))  # shape (n_points, 2)

    return grid_points

    

# Fits grid from make_grid to the adata
def fit_grid_to_adata(grid, dist, kd_tree):
    fit_grid = []
    for point in grid:
        neighbors = kd_tree.query_ball_point(point, dist)
        if len(neighbors) != 0:
            fit_grid.append(point)
    
    return np.array(fit_grid)



## GAUSSIAN KERNEL FUNCTION

# For a given point (x,y), gaussian_kernel_height finds the height of a gaussian kernel 
#   with std. dev. s centered at (0,0)
def gaussian_kernel_height(x, y, s):
    height = np.exp(-((x**2 + y**2) / (2 * s**2)))
    return height



## PARALLELIZED GAUSSIAN APPLICATION FUNCTIONS FOR SMOOTHING

# Apply Gaussian kernel to each point in points_to_smooth (used in both grid-based and in-place smoothing)
def apply_gauss_kernel_chunk(adata_X, 
                             adata_spatial_coords,
                             points_to_smooth,
                             kd_tree,
                             gauss_kernel_function, 
                             min_beads_under_kernel,
                             gaussian_sd,
                             kernel_radius, 
                             total_processes, 
                             process_iter):

    # Figure out what grid_pt range to smooth for current iteration
    pt_start = process_iter * (int(np.ceil(len(points_to_smooth) / total_processes)))
    pt_end = (process_iter + 1) * (int(np.ceil(len(points_to_smooth) / total_processes)))
    
    # Make sure upper range doesn't exceed the number of grid_pts
    if pt_end > len(points_to_smooth):
        pt_end = len(points_to_smooth)

    # Make sure point number is at least 1
    if (pt_end - pt_start) <= 0:
        return -1
    
    # Preallocate the .X array to fill in with smoothed points
    X_sm_chunk = np.empty((pt_end - pt_start, adata_X.shape[1]), dtype=np.float32)

    x_values = adata_spatial_coords[:,0]
    y_values = adata_spatial_coords[:,1]
    
    # For each point in the grid
    for i in range(X_sm_chunk.shape[0]):
        temp_point = points_to_smooth[i + pt_start]
        neighbor_indices = kd_tree.query_ball_point(temp_point, kernel_radius)

        temp_X = adata_X[neighbor_indices, :]
        temp_x_values = x_values[neighbor_indices]
        temp_y_values = y_values[neighbor_indices]

        # Compute the kernel weights
        kernel = gauss_kernel_function(temp_x_values - temp_point[0], 
                                       temp_y_values - temp_point[1], 
                                       gaussian_sd)
        
        # Normalize the kernel
        kernel_sum = np.sum(kernel)
        if kernel_sum > 0:
            kernel /= kernel_sum
        
        # Check if spatial points count under kernel meets threshold
        if len(kernel) < min_beads_under_kernel:
            # Too few points under kernel, set all gene values to 0
            X_sm_chunk[i:i+1, :] = np.zeros(adata_X.shape[1], dtype=np.float32)
        else:
            # Apply the kernel to the gene expression data
            X_sm_chunk[i:i+1, :] = (temp_X.T @ kernel).T.astype(np.float32)

    return X_sm_chunk
    

# Wrapper function compatible with pool.map
def apply_gauss_kernel_wrapper(args):
    return apply_gauss_kernel_chunk(*args)



## GRID-BASED GAUSSIAN SMOOTHING FUNCTION

# Apply grid-based smoothing to ADATA chunk by chunk, 
#   with explicit pool.map parallelization and implicit MM parallelization
def run_gridbased_smoothing(ADATA,
                            STRIDE,
                            GRID_FITTING_DIST,
                            MIN_BEADS_UNDER_KERNEL,
                            GAUSSIAN_SD, 
                            KERNEL_RADIUS,
                            NUM_PROCESSES,
                            NUM_DATA_SPLITS,
                            IS_FFT_EXPT=False):
    t0 = time.time()
    
    # Fit grid points to ADATA
    max_X, min_X = ADATA.obsm['spatial'][:,0].max(), ADATA.obsm['spatial'][:,0].min()
    max_Y, min_Y = ADATA.obsm['spatial'][:,1].max(), ADATA.obsm['spatial'][:,1].min()
    
    # Report ranges of grid
    print("Fitting grid to tissue coordinates...")
    print(f"Variation in X-direction: {max_X - min_X}")
    print(f"Variation in Y-direction: {max_Y - min_Y}")

    # This should always be False, except when hijacking function within select_smoothing_params!
    if IS_FFT_EXPT:
        fit_adata_grid = make_square_grid(max_X, min_X, max_Y, min_Y, STRIDE)
    else:
        # Generate grid
        adata_grid = make_grid(max_X, min_X, max_Y, min_Y, STRIDE)
        
        # Generate KD Tree for efficient grid-fitting
        kd_tree = cKDTree(ADATA.obsm['spatial'])
      
        # Fit grid to adata
        fit_adata_grid = fit_grid_to_adata(adata_grid, GRID_FITTING_DIST, kd_tree)

    # Get x and y ranges of the fit grid 
    x_max_fit_adata_grid = np.max(fit_adata_grid[:,0])
    x_min_fit_adata_grid = np.min(fit_adata_grid[:,0])
    x_range_fit_adata_grid = x_max_fit_adata_grid - x_min_fit_adata_grid
    y_max_fit_adata_grid = np.max(fit_adata_grid[:,1])
    y_min_fit_adata_grid = np.min(fit_adata_grid[:,1])
    y_range_fit_adata_grid = y_max_fit_adata_grid - y_min_fit_adata_grid

    # Report number of grid points
    print(f"Number of points in fitted grid: {len(fit_adata_grid)}")
    
    ## Smooth each spatial chunk of the adata in sequence and combine
    # Initialize lists for batch accumulation
    X_grid_sm_list = []
    fit_adata_grid_list = []
    
    # For each adata chunk
    print(f"Gaussian smoothing will run in {NUM_DATA_SPLITS**2} chunks.")
    print(f"Smoothing chunk: ", end="")
    for i in range(NUM_DATA_SPLITS*NUM_DATA_SPLITS):
        
        # Get current portion of imposed grid points
        x_max_grid_chunk = (
            x_min_fit_adata_grid + (((i % int(NUM_DATA_SPLITS)) + 1) * 
            np.ceil(float(x_range_fit_adata_grid) / float(NUM_DATA_SPLITS)))
        )
        x_min_grid_chunk = (
            x_min_fit_adata_grid + ((i % int(NUM_DATA_SPLITS)) * 
            np.ceil(float(x_range_fit_adata_grid) / float(NUM_DATA_SPLITS)))
        )
        y_max_grid_chunk = (
            y_min_fit_adata_grid + (((i // int(NUM_DATA_SPLITS)) + 1) * 
            np.ceil(float(y_range_fit_adata_grid) / float(NUM_DATA_SPLITS)))
        )
        y_min_grid_chunk = (
            y_min_fit_adata_grid + ((i // int(NUM_DATA_SPLITS)) * 
            np.ceil(float(y_range_fit_adata_grid) / float(NUM_DATA_SPLITS)))
        )

        # Get current portion of adata points (add kernel radius buffer)
        x_max_adata_chunk = x_max_grid_chunk + KERNEL_RADIUS
        x_min_adata_chunk = x_min_grid_chunk - KERNEL_RADIUS
        y_max_adata_chunk = y_max_grid_chunk + KERNEL_RADIUS
        y_min_adata_chunk = y_min_grid_chunk - KERNEL_RADIUS

        mask = ((ADATA.obsm['spatial'][:,0] >= x_min_adata_chunk) &
        (ADATA.obsm['spatial'][:,0] < x_max_adata_chunk) &
        (ADATA.obsm['spatial'][:,1] >= y_min_adata_chunk) &
        (ADATA.obsm['spatial'][:,1] < y_max_adata_chunk))

        # Get sparse chunk of adata.X
        adata_X_chunk = ADATA.X[mask].copy()
        
        # Get spatial coordinates
        adata_spatial_coords_chunk = ADATA.obsm['spatial'][mask]
        
        # Filter fit_adata_grid
        fit_adata_grid_chunk = fit_adata_grid[((fit_adata_grid[:,0] >=  x_min_grid_chunk) &
                                               (fit_adata_grid[:,0] <  x_max_grid_chunk) &
                                               (fit_adata_grid[:,1] >=  y_min_grid_chunk) &
                                               (fit_adata_grid[:,1] <  y_max_grid_chunk))]
        
        # Make new kd_tree
        kd_tree_chunk = cKDTree(adata_spatial_coords_chunk)

        ## Parallelized Smoothing on points within each spatial chunk
        print(f"{i+1}, ", end="")
        with mp.Pool(processes=NUM_PROCESSES) as pool:
            # Run in parallel
            X_sm_results_chunk = pool.map(apply_gauss_kernel_wrapper, 
                                     [(adata_X_chunk, 
                                       adata_spatial_coords_chunk,
                                       fit_adata_grid_chunk,
                                       kd_tree_chunk,
                                       gaussian_kernel_height,
                                       MIN_BEADS_UNDER_KERNEL,
                                       GAUSSIAN_SD, 
                                       KERNEL_RADIUS,
                                       NUM_PROCESSES,
                                       i) for i in range(NUM_PROCESSES)])

        # Free up additional memory
        del adata_X_chunk, adata_spatial_coords_chunk, kd_tree_chunk
        gc.collect()
        time.sleep(0.1)

        # Filter non-empty results and accumulate
        for result in X_sm_results_chunk:
            if isinstance(result, np.ndarray):
                X_grid_sm_list.append(result)
        fit_adata_grid_list.append(fit_adata_grid_chunk)
        
    # After all chunks are done, concatenate all X_grid_sm and fit_adata_grid together
    complete_X_grid_sm = np.concatenate(X_grid_sm_list, axis=0)
    complete_grid_coords = np.concatenate(fit_adata_grid_list, axis=0)

    del X_grid_sm_list, fit_adata_grid_list
    gc.collect()
    time.sleep(0.1)

    # Make final sm_adata
    complete_sm_adata = ad.AnnData(complete_X_grid_sm.astype(np.float32))
    complete_sm_adata.obsm['spatial'] = complete_grid_coords

    del complete_X_grid_sm, complete_grid_coords
    gc.collect()
    time.sleep(0.1)

    # Add .var(m) gene info and general .uns info to final sm_adata
    complete_sm_adata.var = ADATA.var
    complete_sm_adata.varm = ADATA.varm
    complete_sm_adata.uns = ADATA.uns

    runtime = time.time() - t0
    print("\nTotal runtime for grid-based Gaussian smoothing: {:.2f} seconds.\n".format(runtime))
                                                         
    return complete_sm_adata



## IN-PLACE GAUSSIAN SMOOTHING FUNCTION
    
# Apply in-place smoothing to ADATA chunk by chunk, 
#   with explicit pool.map parallelization and implicit MM parallelization
def run_inplace_smoothing(ADATA,
                          MIN_BEADS_UNDER_KERNEL,
                          GAUSSIAN_SD, 
                          KERNEL_RADIUS,
                          NUM_PROCESSES,
                          NUM_DATA_SPLITS):
    t0 = time.time()
    
    # Get x and y ranges of the fit grid 
    max_x, min_x = ADATA.obsm['spatial'][:,0].max(), ADATA.obsm['spatial'][:,0].min()
    max_y, min_y = ADATA.obsm['spatial'][:,1].max(), ADATA.obsm['spatial'][:,1].min()
    range_x, range_y = (max_x - min_x), (max_y - min_y) 
    
    # initialize lists to store smoothed adata.X chunks and adata.obs and adata.obsm chunks
    X_sm_list = []
    obs_list = []
    obsm_list = []
    
    # For each chunk of the adata
    print(f"Gaussian smoothing will run in {NUM_DATA_SPLITS**2} chunks.")
    print(f"Smoothing chunk: ", end="")
    for i in range(NUM_DATA_SPLITS*NUM_DATA_SPLITS):
        
        # Get chunk of adata points
        x_max_chunk = (
            min_x + (((i % int(NUM_DATA_SPLITS)) + 1) * 
            np.ceil(float(range_x) / float(NUM_DATA_SPLITS)))
        )
        x_min_chunk = (
            min_x + ((i % int(NUM_DATA_SPLITS)) * 
            np.ceil(float(range_x) / float(NUM_DATA_SPLITS)))
        )
        y_max_chunk = (
            min_y + (((i // int(NUM_DATA_SPLITS)) + 1) * 
            np.ceil(float(range_y) / float(NUM_DATA_SPLITS)))
        )
        y_min_chunk = (
            min_y + ((i // int(NUM_DATA_SPLITS)) * 
            np.ceil(float(range_y) / float(NUM_DATA_SPLITS)))
        )

        # Get chunk of adata points with KERNEL_RADIUS padding
        x_max_padded_chunk = x_max_chunk + KERNEL_RADIUS
        x_min_padded_chunk = x_min_chunk - KERNEL_RADIUS
        y_max_padded_chunk = y_max_chunk + KERNEL_RADIUS
        y_min_padded_chunk = y_min_chunk - KERNEL_RADIUS

        yes_padding_mask = ((ADATA.obsm['spatial'][:,0] >= x_min_padded_chunk) &
                            (ADATA.obsm['spatial'][:,0] < x_max_padded_chunk) &
                            (ADATA.obsm['spatial'][:,1] >= y_min_padded_chunk) &
                            (ADATA.obsm['spatial'][:,1] < y_max_padded_chunk))

        no_padding_mask = ((ADATA.obsm['spatial'][:,0] >= x_min_chunk) &
                           (ADATA.obsm['spatial'][:,0] < x_max_chunk) &
                           (ADATA.obsm['spatial'][:,1] >= y_min_chunk) &
                           (ADATA.obsm['spatial'][:,1] < y_max_chunk))

        # Filter sparse adata.X for chunk with KERNEL_RADIUS padding
        adata_X_chunk = ADATA.X[yes_padding_mask].copy()

        # Filter chunk coordinates with KERNEL_RADIUS padding
        adata_padded_coords_chunk = ADATA.obsm['spatial'][yes_padding_mask]
        
        # Filter chunk coordinates with no padding
        adata_coords_chunk = ADATA.obsm['spatial'][no_padding_mask]
        
        # Make new kd_tree
        kd_tree_chunk = cKDTree(adata_padded_coords_chunk)

        ## Parallelized Smoothing on points within each spatial chunk
        print(f"{i+1}, ", end="")
        with mp.Pool(processes=NUM_PROCESSES) as pool:
            X_sm_results_chunk = pool.map(apply_gauss_kernel_wrapper, 
                                         [(adata_X_chunk, 
                                           adata_padded_coords_chunk,
                                           adata_coords_chunk,
                                           kd_tree_chunk,
                                           gaussian_kernel_height,
                                           MIN_BEADS_UNDER_KERNEL,
                                           GAUSSIAN_SD, 
                                           KERNEL_RADIUS,
                                           NUM_PROCESSES,
                                           i) for i in range(0,NUM_PROCESSES)])

        # Free up chunk memory
        del adata_X_chunk, adata_coords_chunk, kd_tree_chunk
        gc.collect()
        time.sleep(0.1)

        # Filter non-empty results and accumulate
        for result in X_sm_results_chunk:
            if isinstance(result, np.ndarray):
                X_sm_list.append(result)

        # Filter chunk .obs and .obsm for reconstructing final smoothed adata
        adata_obs_chunk = ADATA.obs[no_padding_mask]
        adata_obsm_chunk = {key: value[no_padding_mask] for key, value in ADATA.obsm.items()}  
        obs_list.append(adata_obs_chunk)
        obsm_list.append(adata_obsm_chunk)
        
    # After all chunks are done, concatenate all .X, .obs, and .obsm lists together
    complete_X_sm = np.concatenate(X_sm_list, axis=0)
    complete_obs = pd.concat(obs_list, axis=0)
    complete_obsm = {key: np.concatenate([obsm[key] for obsm in obsm_list], axis=0) for key in obsm_list[0]}

    del X_sm_list, obs_list, obsm_list
    gc.collect()
    time.sleep(0.1)

    # Make final sm_adata
    complete_sm_adata = ad.AnnData(complete_X_sm.astype(np.float32))
    
    complete_sm_adata.obs = complete_obs
    complete_sm_adata.obsm = complete_obsm
    complete_sm_adata.var = ADATA.var
    complete_sm_adata.varm = ADATA.varm
    complete_sm_adata.uns = ADATA.uns

    del complete_X_sm, complete_obs, complete_obsm
    gc.collect()
    time.sleep(0.1)

    runtime = time.time() - t0
    print("\nTotal runtime for in-place Gaussian smoothing: {:.2f} seconds.\n".format(runtime))
                                                         
    return complete_sm_adata


def perform_density_check(adata, grid_based_or_not, gaussian_sd, min_spots_under_gaussian, 
                          stride=None, grid_fitting_dist=None, num_data_splits=None):
    """
    Helper function to pre-scan data in spatial chunks to estimate drop-out rates.
    Uses chunking + padding to keep memory low, mirroring the smoothing logic.
    """
    print(f"Checking dataset point density...")

    # Define layout parameters
    kernel_radius = 3 * gaussian_sd
    max_x, min_x = adata.obsm['spatial'][:,0].max(), adata.obsm['spatial'][:,0].min()
    max_y, min_y = adata.obsm['spatial'][:,1].max(), adata.obsm['spatial'][:,1].min()
    range_x, range_y = (max_x - min_x), (max_y - min_y)
    
    # Initialize counters
    total_ignored = 0
    total_checked = 0
    
    # If splits not provided, estimate based on size (same logic as main func)
    if num_data_splits is None:
        locs = adata.shape[0]
        if float(locs) < 100000: num_data_splits = 2
        elif float(locs) < 1000000: num_data_splits = 4
        else: num_data_splits = 6

    # Iterate through chunks
    for i in range(num_data_splits * num_data_splits):
        
        # 1. Define Chunk Boundaries (Unpadded) - The region we want to TEST
        x_min_chunk = min_x + ((i % num_data_splits) * np.ceil(range_x / num_data_splits))
        x_max_chunk = min_x + (((i % num_data_splits) + 1) * np.ceil(range_x / num_data_splits))
        y_min_chunk = min_y + ((i // num_data_splits) * np.ceil(range_y / num_data_splits))
        y_max_chunk = min_y + (((i // num_data_splits) + 1) * np.ceil(range_y / num_data_splits))

        # 2. Define Padded Boundaries - The region we build the TREE on (Reference)
        # We need padding so points at the edge of the chunk can find their neighbors
        x_min_pad = x_min_chunk - kernel_radius
        x_max_pad = x_max_chunk + kernel_radius
        y_min_pad = y_min_chunk - kernel_radius
        y_max_pad = y_max_chunk + kernel_radius

        # 3. Filter Data for Reference Tree (Padded)
        # These are the neighbors we might find
        mask_pad = ((adata.obsm['spatial'][:,0] >= x_min_pad) & (adata.obsm['spatial'][:,0] < x_max_pad) &
                    (adata.obsm['spatial'][:,1] >= y_min_pad) & (adata.obsm['spatial'][:,1] < y_max_pad))
        
        coords_pad = adata.obsm['spatial'][mask_pad]
        
        # If chunk is empty, skip
        if len(coords_pad) == 0:
            continue

        # 4. Build Local Tree (Small & Fast)
        chunk_tree = cKDTree(coords_pad)

        # 5. Define Query Points (Unpadded)
        # These are the points we are actually testing for density
        if grid_based_or_not:
            # Generate grid for this chunk only
            chunk_grid = make_grid(x_max_chunk, x_min_chunk, y_max_chunk, y_min_chunk, stride)
            if len(chunk_grid) == 0: continue
            
            # Filter grid points: Must be within grid_fitting_dist of ACTUAL data
            # We reuse the chunk_tree for this!
            valid_grid_points = fit_grid_to_adata(chunk_grid, grid_fitting_dist, chunk_tree)
            query_points = valid_grid_points
            
        else:
            # For in-place, query points are the actual data points in the unpadded region
            mask_no_pad = ((adata.obsm['spatial'][:,0] >= x_min_chunk) & (adata.obsm['spatial'][:,0] < x_max_chunk) &
                           (adata.obsm['spatial'][:,1] >= y_min_chunk) & (adata.obsm['spatial'][:,1] < y_max_chunk))
            query_points = adata.obsm['spatial'][mask_no_pad]

        if len(query_points) == 0:
            continue

        # 6. Subsample Query Points
        # We don't need to check every point to get a percentage. 
        # Capping at 100k per chunk is very safe and fast.
        MAX_SAMPLES_PER_CHUNK = 100_000
        if len(query_points) > MAX_SAMPLES_PER_CHUNK:
            indices = np.random.choice(len(query_points), MAX_SAMPLES_PER_CHUNK, replace=False)
            query_subset = query_points[indices]
        else:
            query_subset = query_points

        # 7. Query Tree for Density
        # How many neighbors does each query point have within kernel_radius?
        try:
            # return_length=True is much faster (requires Scipy >= 1.8)
            neighbor_counts = chunk_tree.query_ball_point(query_subset, kernel_radius, return_length=True)
        except TypeError:
            # Fallback
            neighbor_lists = chunk_tree.query_ball_point(query_subset, kernel_radius)
            neighbor_counts = np.array([len(l) for l in neighbor_lists])

        if not isinstance(neighbor_counts, np.ndarray):
            neighbor_counts = np.array(neighbor_counts)

        # 8. Update Tally
        ignored_in_chunk = np.sum(neighbor_counts < min_spots_under_gaussian)
        
        total_ignored += ignored_in_chunk
        total_checked += len(query_subset)

        # delete large variables
        del coords_pad, chunk_tree, query_points, query_subset, neighbor_counts

    # 9. Report Findings
    ignored_percent = (total_ignored / total_checked) * 100
    
    if ignored_percent > 15:
        raise ValueError(
            f"Error: Approx {ignored_percent:.2f}% of points "
            f"will be ignored due to min_spots_under_gaussian={min_spots_under_gaussian}. \n"
            "A large fraction of data is ignored due to the choice of gaussian_sd and min_spots_under_gaussian. "
            "Please double check that data has gone through quality control to filter low UMI spots. Please double check your choice of gaussian_sd. You may fix this error by choosing a lower value for min_spots_under_gaussian."
        )
    elif ignored_percent > 2:
        print(f"Warning: Approx {ignored_percent:.2f}% of points will be discarded due to low neighbors.")
        

## EXECUTABLE FUNCTION

def run_parallelized_smoothing(adata,
                               grid_based_or_not,
                               gaussian_sd,
                               min_spots_under_gaussian,
                               stride=None,
                               grid_fitting_dist=None,
                               num_processes=10,
                               num_data_splits=None):
    """
    Performs parallelized Gaussian smoothing on spatial transcriptomics data.

    Parameters:
    -----------
    adata : AnnData
        The input dataset (an AnnData object).
    grid_based_or_not : bool
        True = grid-based smoothing, smooth only at imposed hexagonal grid points.
        False = in-place smoothing, smooth at every spatial location in the dataset.
        In-place is effective for data where spatial spots have "cell-sized" or larger resolution (10-50 micron).
        Grid-based is effective for data where spatial spots have subcellular to "cell-sized" resolution (0.5-10 micron).
    gaussian_sd : float
        Standard deviation for the Gaussian kernel. !! Carefully choose this variable based on S.T. data platform. !!
        * gaussian_sd = target_microns * micron_to_unit_conversion *
        - target_microns = 20 to 40 is generally ideal across high-resolution S.T. platforms
        - micron_to_unit_conversion is different for each S.T. data platform
        --> Check Smoothie/docs/guides/micron_to_unit_conversion_table.md for these values.
    min_spots_under_gaussian : int
        Minimum number of data points within radius (3 * gaussian_sd) of center point 
            for smoothing to occur at that location (Default is 25-100).
        If chosen min_spots_under_gaussian is too high for the coordinate density, a warning message will appear.
    stride : float, optional
        Stride value for grid-based smoothing (default stride = 1 * gaussian_sd). 
            (0.5 * gaussian_sd allows for a denser grid).
    grid_fitting_dist : float
        The minimum distance a grid point must be from a spatial coordinate to be retained in the grid.
        (Default is 0.25 * gaussian_sd. Values up to 0.5 * gaussian_sd are reasonable too).
    num_processes : int, optional
        Number of parallel processes to use (Default is 10).
    num_data_splits : int, optional
        Data will be smoothed in num_data_splits x num_data_splots chunks for memory efficiency. 
        If unspecified, it will automatically be selected for speed and memory balance.

    Returns:
    --------
    sm_adata : AnnData
        The smoothed AnnData object.
    """

    ## ASSERTIONS FOR CODE

    # adata
    assert hasattr(adata, "X"), "adata.X does not exist"
    assert adata.X.ndim == 2, "adata.X must be a 2D matrix"
    assert 'spatial' in adata.obsm, "'spatial' key is missing in adata.obsm"
    adata.obsm['spatial'] = np.array(adata.obsm['spatial'])
    assert np.array(adata.obsm['spatial']).ndim == 2, "adata.obsm['spatial'] must be a 2D array"

    # Check for unique gene names
    if not adata.var_names.is_unique:
        raise ValueError(
            "Gene names (adata.var_names) are not unique. "
            "Please run `adata.var_names_make_unique()` before running this function."
        )

    # grid_based_or_not
    assert isinstance(grid_based_or_not, bool), "grid_based_or_not must be True or False"

    # gaussian_sd
    assert isinstance(gaussian_sd, (int, float)) and gaussian_sd > 0, "gaussian_sd must be a positive number."

    # min_spots_under_gaussian
    assert isinstance(min_spots_under_gaussian, int) and min_spots_under_gaussian > 0, \
    "min_spots_under_gaussian must be a positive integer. (Range ~25-100)."

    # stride
    if grid_based_or_not:
        if stride is None:
            stride = gaussian_sd
        assert isinstance(stride, (int, float)) and stride >= (0.5*gaussian_sd) and stride <= gaussian_sd, (
            "stride must be a value in the range [0.5 * gaussian_sd, 1.0 * gaussian_sd]. (Default is 1.0 * gaussian_sd)."
        )

    # grid_fitting_dist
    if grid_based_or_not:
        if grid_fitting_dist is None:
            grid_fitting_dist = 0.25 * gaussian_sd
        assert (isinstance(grid_fitting_dist, (int, float))
                and grid_fitting_dist >= (0.25*gaussian_sd) 
                and grid_fitting_dist <= (0.5*gaussian_sd)), (
            "grid_fitting_dist must be a value in the range [0.25 * gaussian_sd, 0.5 * gaussian_sd]. (Default is 0.25 * gaussian_sd)."
        )
               
    # num_processes
    assert isinstance(num_processes, int) and num_processes > 0, (
        "num_processes must be a positive integer. (Default is 10)."
    )

    # num_data_splits
    if num_data_splits is not None:
        assert isinstance(num_data_splits, int) and num_data_splits > 0, (
            "num_data_splits must be a positive integer. (Default is 1 to 6)."
        )
    if num_data_splits is None:
        locs = adata.shape[0]  # Get number of locations
        if float(locs) < 100000:
            num_data_splits = 2
        elif float(locs) < 1000000:
            num_data_splits = 4
        else:
            num_data_splits = 6

    
    ## RUN SMOOTHING 
    
    # ensure adata.X is a sparse csr matrix
    if issparse(adata.X):
        adata.X = adata.X.tocsr()
    else:
        adata.X = csr_matrix(adata.X)

    # Set default value of kernel_radius
    kernel_radius = 3 * gaussian_sd
    
    # If grid-based smoothing
    if grid_based_or_not:

        # Check value for min_spots_under_gaussian
        perform_density_check(adata, 
                              grid_based_or_not, 
                              gaussian_sd, 
                              min_spots_under_gaussian,
                              stride,
                              grid_fitting_dist, 
                              num_data_splits)
        
        sm_adata = run_gridbased_smoothing(adata,
                                           stride,
                                           grid_fitting_dist,
                                           min_spots_under_gaussian,
                                           gaussian_sd, 
                                           kernel_radius,
                                           num_processes,
                                           num_data_splits)
    # If in-place smoothing
    else:
        locs, genes = adata.shape  # Get number of locations and genes
        memory_estimation_gb = float(locs * genes * 4) / 1e9
        print(f"Storage for smoothed count matrix will be {memory_estimation_gb:.2f} GB.")
        
        # Assert that memory estimate is manageable
        assert memory_estimation_gb < 100, "Memory usage will be high (>100GB)! Consider grid-based smoothing."

        # Check value for min_spots_under_gaussian
        perform_density_check(adata, 
                              grid_based_or_not, 
                              gaussian_sd, 
                              min_spots_under_gaussian,
                              stride=None,
                              grid_fitting_dist=None, 
                              num_data_splits=num_data_splits)
        
        sm_adata = run_inplace_smoothing(adata,
                                         min_spots_under_gaussian,
                                         gaussian_sd, 
                                         kernel_radius,
                                         num_processes,
                                         num_data_splits)

    # Post-smoothing cleanup
    gene_sums = np.sum(sm_adata.X, axis=0)
    zero_genes_mask = gene_sums == 0
    
    if np.any(zero_genes_mask):
        dropped_genes = sm_adata.var_names[zero_genes_mask]
        print("The following genes were removed because of weak signal in the dataset:")
        for gene in dropped_genes:
            print(f"{gene}")
        sm_adata = sm_adata[:, ~zero_genes_mask]
        
    ## RETURN SMOOTHED ADATA
    return sm_adata