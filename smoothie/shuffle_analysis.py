# Copyright (c) 2026 Chase Holdener
# Licensed under the MIT License. See LICENSE file for details.

"""
Memory-efficient spatial shuffling and correlation analysis.

This module provides functions for null distribution generation through spatial shuffling,
with careful memory management to handle large datasets.
"""

import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix
import gc
from typing import Union, List, Tuple


def compute_shuffled_correlation_percentiles(
    adata: Union[ad.AnnData, List[ad.AnnData]],
    grid_based_or_not,
    gaussian_sd,
    min_spots_under_gaussian,
    stride: float = None,
    grid_fitting_dist: float = None,
    num_processes: int = 4,
    num_data_splits: int = None,
    seed: int = None
) -> Tuple[float, float, float]:
    """
    Shuffles spatial coordinates, smoothes data, and computes correlation matrix percentiles.
    
    Parameters
    ----------
    adata : ad.AnnData or List[ad.AnnData]
        Single AnnData object or list of AnnData objects to process.
        For multiple datasets, genes will be aligned.
    grid_based_or_not : bool
        True = bin shuffling + grid-based smoothing (for subcellular resolution data)
        False = in-place shuffling + in-place smoothing (for cellular resolution data)
    gaussian_sd : float
        gaussian_sd = target_microns * micron_to_unit_conversion value
        Standard deviation of Gaussian kernel for smoothing (in coordinate units).
    min_spots_under_gaussian : int
        Minimum number of spots required under kernel for valid smoothing.
    stride : float, optional
        Grid spacing for smoothing (default = gaussian_sd).
        Only used if grid_based_or_not=True.
    grid_fitting_dist : float, optional
        Minimum distance from grid point to data for grid fitting
        (default = 0.25 × gaussian_sd). Only used if grid_based_or_not=True.
    num_processes : int, default=4
        Number of parallel processes for smoothing.
    num_data_splits : int, optional
        Number of data splits for parallel processing.
        If None, automatically determined for memory efficiency.
    seed : int, optional
        Random seed for reproducibility.
    
    Returns
    -------
    p95 : float
        95th percentile of correlation coefficients (lower triangle).
    p99 : float
        99th percentile of correlation coefficients (lower triangle).
    p999 : float
        99.9th percentile of correlation coefficients (lower triangle).
    
    Notes
    -----
    - Call this function with the SAME parameters that were used for run_parallelized_smoothing.
    """
    # Import smoothie functions (avoid circular import)
    from .gaussian_smoothing import bin_shuffle_adata, run_parallelized_smoothing
    from . import compute_correlation_matrix
    
    # Set defaults
    if stride is None:
        stride = gaussian_sd
    if grid_fitting_dist is None:
        grid_fitting_dist = 0.25 * gaussian_sd
    kernel_radius = 3 * gaussian_sd
    
    if seed is not None:
        np.random.seed(seed)
    
    # Handle single AnnData vs list
    is_single = not isinstance(adata, list)
    adata_list = [adata] if is_single else adata
    
    # Derive shuffle method from grid_based_or_not
    shuffle_method = 'bin' if grid_based_or_not else 'inplace'
    
    print(f"Running shuffled data spatial correlation testing...")
    
    # Step 1: Shuffle coordinates for each dataset (in-place to save memory)
    print("\n" + "=" * 60 + f"\n[Step 1/4] Shuffling coordinates...")
    shuffled_adata_list = []
    
    for i, adata_orig in enumerate(adata_list):
        print(f"Dataset {i+1}/{len(adata_list)}, ", end="")
        
        # Create minimal copy with only essential data
        adata_shuffled = ad.AnnData(
            X=adata_orig.X.copy() if hasattr(adata_orig.X, 'copy') else adata_orig.X,
            var=adata_orig.var[[]].copy(),  # Empty DataFrame with same index
            obs=adata_orig.obs[[]].copy()   # Empty DataFrame with same index
        )
        adata_shuffled.var_names = adata_orig.var_names.copy()
        adata_shuffled.obsm['spatial'] = adata_orig.obsm['spatial'].copy()
        
        if shuffle_method == 'inplace':
            # Simple random permutation of all coordinates
            np.random.shuffle(adata_shuffled.obsm['spatial'])
            
        elif shuffle_method == 'bin':
            # Shuffle spatial bins of points (maintaining structure within each bin)
            bin_shuffle_adata(adata_shuffled, stride, seed)
           
        else:
            raise ValueError(f"grid_based_or_not must be True or False.")
        
        shuffled_adata_list.append(adata_shuffled)
    
    # Step 2: Smooth each shuffled dataset
    print("\n\n" + "=" * 60 + f"\n[Step 2/4] Smoothing shuffled datasets...")
    smoothed_adata_list = []
    
    for i, adata_shuffled in enumerate(shuffled_adata_list):
        print(f"Dataset {i+1}/{len(shuffled_adata_list)}:")
        
        sm_sh_adata = run_parallelized_smoothing(
            adata=adata_shuffled,
            grid_based_or_not=grid_based_or_not,
            gaussian_sd=gaussian_sd,
            min_spots_under_gaussian=min_spots_under_gaussian,
            stride=stride if grid_based_or_not else None,
            grid_fitting_dist=grid_fitting_dist if grid_based_or_not else None,
            num_processes=num_processes,
            num_data_splits=num_data_splits
        )
        
        # Keep only essential data
        sm_sh_adata_minimal = ad.AnnData(
            X=sm_sh_adata.X,
            var=sm_sh_adata.var[[]].copy()
        )
        sm_sh_adata_minimal.var_names = sm_sh_adata.var_names.copy()
        
        smoothed_adata_list.append(sm_sh_adata_minimal)
        
        # Clean up
        del sm_sh_adata
        del adata_shuffled
        gc.collect()
    
    # Delete original shuffled data
    del shuffled_adata_list
    gc.collect()
    
    # Step 3: Align genes if multiple datasets
    if len(smoothed_adata_list) > 1:
        print("=" * 60 + f"\n[Step 3/4] Aligning genes across datasets...")
        
        # Find UNION of genes
        all_genes = set()
        for sm_sh_adata in smoothed_adata_list:
            all_genes.update(sm_sh_adata.var_names)
        all_genes_list = sorted(list(all_genes))
        
        print(f"Total unique genes across datasets: {len(all_genes_list)}")
        
        # Create mapping: gene_name -> global_index
        gene_to_idx = {gene: i for i, gene in enumerate(all_genes_list)}
        
        # Process datasets sequentially and concatenate on-the-fly
        X_combined = None
        
        print("Merging datasets: ", end="")
        
        for i, sm_sh_adata in enumerate(smoothed_adata_list):
            print(f"{i+1}, ", end="", flush=True)
            
            # 1. Get data as dense numpy array
            curr_X = sm_sh_adata.X
            
            # 2. Create Aligned Matrix
            # Initialize with zeros (implicitly handles missing genes)
            n_cells = curr_X.shape[0]
            n_total_genes = len(all_genes_list)
            
            # Use same dtype as input (usually float32) to save memory
            X_aligned = np.zeros((n_cells, n_total_genes), dtype=curr_X.dtype)
            
            # 3. Map local columns to global columns
            local_genes = sm_sh_adata.var_names
            global_indices = [gene_to_idx[g] for g in local_genes]
            
            # 4. Fill values
            X_aligned[:, global_indices] = curr_X
            
            # 5. Concatenate
            if X_combined is None:
                X_combined = X_aligned
            else:
                X_combined = np.concatenate([X_combined, X_aligned], axis=0)
                
            # Clean up immediately
            del X_aligned, curr_X
            gc.collect()
        
        # Clean up smoothed_adata_list
        del smoothed_adata_list
        gc.collect()
        
        # Create final combined AnnData
        combined_adata = ad.AnnData(X=X_combined)
        combined_adata.var_names = all_genes_list
        
        del X_combined
        gc.collect()
        
        final_adata = combined_adata
        
    else:
        print("\n" + "=" * 60 + f"\n[Step 3/4] Single dataset - no gene alignment needed")
        final_adata = smoothed_adata_list[0]
             
        del smoothed_adata_list
        gc.collect()
    
    # Step 4: Compute correlation matrix and extract percentiles
    print("\n\n" + "=" * 60 + f"\n[Step 4/4] Computing correlation matrix and percentiles...")
    print(f"Matrix size: {final_adata.n_vars} × {final_adata.n_vars} genes")
    
    # Compute correlation matrix
    pearsonR_mat, _ = compute_correlation_matrix(final_adata.X)
    
    # Clean up AnnData - no longer needed
    del final_adata
    gc.collect()
    
    # Extract lower triangle (excluding diagonal)
    n = pearsonR_mat.shape[0]
    lower_tri_indices = np.tril_indices(n, k=-1)
    lower_tri_values = pearsonR_mat[lower_tri_indices]
    
    # Compute percentiles
    p95 = np.percentile(lower_tri_values, 95)
    p99 = np.percentile(lower_tri_values, 99)
    p999 = np.percentile(lower_tri_values, 99.9)
    
    # Clean up
    del pearsonR_mat, lower_tri_values
    gc.collect()
    
    # Report results
    print("\n" + "=" * 60 + f"\nResults: Shuffled Correlation Percentiles")
    print(f"95th percentile:  {p95:.6f}")
    print(f"99th percentile:  {p99:.6f}")
    print(f"99.9th percentile: {p999:.6f}")
    
    return p95, p99, p999