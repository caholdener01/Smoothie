# Copyright (c) 2025 Chase Holdener
# Licensed under the MIT License. See LICENSE file for details.

import numpy as np
import pandas as pd
import anndata as ad
from scipy.stats import pearsonr


def concatenate_smoothed_matrices(sm_adata_set):
    """
    Concatenates smoothed count matrices from multiple datasets, aligning gene columns.
    
    Parameters:
    sm_adata_set (list of AnnData): List of AnnData objects with smoothed count matrices.
    
    Returns:
    tuple: 
        - X_concat (np.ndarray): Concatenated matrix with genes aligned across datasets.
        - gene_names (list of str): List of gene names, each corresponding to a col of X_concat.
    """
    
    # Collect unique gene names and sort them
    gene_names = sorted(set.union(*(set(sm_adata.var_names) for sm_adata in sm_adata_set)))
    gene_names_to_index = {gene: idx for idx, gene in enumerate(gene_names)}

    # Determine total obs and precompute obs index ranges
    obs_counts = np.array([sm_adata.shape[0] for sm_adata in sm_adata_set])
    sm_adata_obs_ranges = np.concatenate(([0], np.cumsum(obs_counts)))

    # Preallocate matrix
    X_concat = np.zeros((sm_adata_obs_ranges[-1], len(gene_names)), dtype=np.float32)

    # Fill the matrix
    for i, sm_adata in enumerate(sm_adata_set):
        row_start, row_end = sm_adata_obs_ranges[i], sm_adata_obs_ranges[i+1]
        col_indices = [gene_names_to_index[gene] for gene in sm_adata.var_names]
        X_concat[row_start:row_end, col_indices] = sm_adata.X

    return X_concat, gene_names


def align_gene_correlation_matrices(sm_adata_set, pearsonR_mat_set):
    """
    Aligns gene correlation matrices such that each row and column index corresponds to the same gene across all matrices.
    
    Parameters:
    sm_adata_set (list of AnnData): List of AnnData objects with smoothed count matrices.
    pearsonR_mat_set (list of np.ndarray): List of Pearson correlation matrices corresponding to each dataset.
    
    Returns:
    aligned_pearsonR_tensor (np.ndarray): A 3D tensor where each slice along the third axis is an aligned Pearson correlation matrix.
    """
    
    # Collect unique gene names and sort them
    gene_names = sorted(set.union(*(set(sm_adata.var_names) for sm_adata in sm_adata_set)))
    gene_names_to_index = {gene: idx for idx, gene in enumerate(gene_names)}
    
    # Initialize a 3D tensor of pearsonR matrices for each dataset as np.nan
    aligned_pearsonR_tensor = np.full((len(gene_names), len(gene_names), len(sm_adata_set)), 
                                      np.nan, dtype=np.float32)
    
    # Fill in the tensor for each dataset
    for d, sm_adata in enumerate(sm_adata_set):
        # Get indices of genes that exist in gene_names
        mask = np.array([gene in gene_names_to_index for gene in sm_adata.var_names])
        valid_genes = np.array(sm_adata.var_names)[mask]
        
        # Get corresponding indices
        idx_in_global = np.array([gene_names_to_index[gene] for gene in valid_genes])
        idx_in_local = np.where(mask)[0]

        # Extract the submatrix from the current dataset
        submatrix = pearsonR_mat_set[d][idx_in_local[:, None], idx_in_local]
        
        # Assign values using NumPy advanced indexing
        aligned_pearsonR_tensor[idx_in_global[:, None], idx_in_global, d] = submatrix 

    return aligned_pearsonR_tensor


def create_gene_embeddings(aligned_pearsonR_tensor, sm_adata_set, pcc_cutoff, node_label_df=None, E_max=25, seed=0):
    """
    Extracts gene embeddings by filtering aligned Pearson correlation matrices for spatially informative genes  
    and optionally balancing embedding features based on gene module assignments.
    
    Parameters:
    -----------
    aligned_pearsonR_tensor (np.ndarray):  
        A 3D tensor containing aligned Pearson correlation matrices across datasets.  
    
    sm_adata_set (list of AnnData):  
        List of AnnData objects containing smoothed count matrices.  
    
    pcc_cutoff (float):  
        Minimum Pearson correlation coefficient a gene must exceed in at least one dataset  
        to be included in the embedding vectors.  
    
    node_label_df (pd.DataFrame, optional):  
        DataFrame containing gene module assignments. If provided, gene features will be balanced  
        across modules to ensure diversity in embedding features. Default is None.  
    
    E_max (int, optional):  
        Maximum number of genes per module that can be used as embedding features.  
        Only applies if `node_label_df` is provided. Default is 25.  
    
    seed (int, optional):  
        Random seed for reproducibility of embedding feature downsampling. Default is 0.  
    
    Returns:
    --------
    tuple:
        gene_embeddings_tensor (np.ndarray):  
            A 3D tensor of gene embeddings with the following dimensions:  
            - Axis 0: Gene embedding vectors for comparison.  
            - Axis 1: Embedding features comprising the vectors.  
            - Axis 2: Spatial datasets, in the same order as `sm_adata_set`.  
        
        robust_gene_names (np.ndarray):  
            Array of gene names corresponding to the first axis (embedding vectors) of `gene_embeddings_tensor`.  
    """
    
    # Collect unique gene names and sort them
    gene_names = sorted(set.union(*(set(sm_adata.var_names) for sm_adata in sm_adata_set)))
    
    # Compute max correlation for each gene across datasets
    for d in range(aligned_pearsonR_tensor.shape[2]):
        np.fill_diagonal(aligned_pearsonR_tensor[:, :, d], 0.0)
        
    gene_max = np.nanmax(aligned_pearsonR_tensor, axis=(1, 2))

    for d in range(aligned_pearsonR_tensor.shape[2]):
        np.fill_diagonal(aligned_pearsonR_tensor[:, :, d], 1.0)

    # Filter genes above the PCC cutoff
    gene_above_cutoff = gene_max > pcc_cutoff
    robust_gene_names = np.array(gene_names)[gene_above_cutoff]

    # Extract submatrix for robust genes
    aligned_pearsonR_tensor_robust = aligned_pearsonR_tensor[gene_above_cutoff, :, :][:, gene_above_cutoff, :]
    epsilon = 1e-5  # Small value to avoid division by zero
    aligned_pearsonR_tensor_robust = np.clip(aligned_pearsonR_tensor_robust, -1 + epsilon, 1 - epsilon)
    
    # Convert to Fisher’s Z-Scores
    Zscore_aligned_pearsonR_tensor_robust = 0.5 * np.log((1 + aligned_pearsonR_tensor_robust) / (1 - aligned_pearsonR_tensor_robust))

    # Downsample large gene module embedding features if modules are provided
    if node_label_df is not None and E_max is not None:
        def sample_indices(group):
            return np.random.choice(group.index, min(len(group), E_max), replace=False)

        np.random.seed(seed)
        sampled_indices = node_label_df.groupby('community_label', group_keys=False).apply(sample_indices).explode().astype(int)
        
        # Create a mask for selected genes
        selected_genes = set(node_label_df.loc[sampled_indices, 'name'])
        downsampling_mask = np.array([gene in selected_genes for gene in robust_gene_names])

        # Create gene embeddings with downsampling
        gene_embeddings_tensor = Zscore_aligned_pearsonR_tensor_robust[:, downsampling_mask, :]

    else: # Take gene embeddings without module-based feature downsampling 
        gene_embeddings_tensor = Zscore_aligned_pearsonR_tensor_robust

    # Return gene embeddings and the gene names for each row vector
    return gene_embeddings_tensor, robust_gene_names


def get_gene_stabilities_across_datasets(gene_embeddings_tensor, robust_gene_names, sm_adata_set):
    """
    Computes stability of same gene embeddings across dataset pairs using Pearson correlation.
    
    Parameters:
    -----------
    gene_embeddings_tensor (np.ndarray): 
        Tensor of gene embeddings output from function create_gene_embeddings.
    
    robust_gene_names (np.ndarray):
        Array of gene names corresponding to axis 0 of `gene_embeddings_tensor`.
    
    sm_adata_set (list of AnnData):
        List of AnnData objects containing smoothed count matrices.
    
    Returns:
    --------
    gene_stabilities_across_datasets (np.ndarray):
        A tensor containing stability values for each gene across dataset pairs.
        - Axis 0: Spatial datasets, in the same order as `sm_adata_set`
        - Axis 1: Spatial datasets, in the same order as `sm_adata_set`
        - Axis 2: Genes, in same order as robust_gene_names.

        gene_stabilities_across_datasets[i,j,k] is the stability of gene
        robust_gene_names[k] between dataset sm_adata_set[i] and dataset sm_adata_set[j] !
    """
    
    num_genes, num_features, num_datasets = gene_embeddings_tensor.shape
    
    # Initialize results tensor with NaNs
    gene_stabilities_across_datasets = np.full((num_datasets, num_datasets, num_genes), np.nan, dtype=np.float32)

    # Iterate over genes
    for g in range(num_genes):
        goi_embeddings = gene_embeddings_tensor[g]

        # Compute valid (non-NaN, non-Inf) masks
        valid_masks = ~(np.isnan(goi_embeddings) | np.isinf(goi_embeddings))

        # Compute pairwise correlations
        for i in range(num_datasets):
            for j in range(i):
                mask = valid_masks[:, i] & valid_masks[:, j]
                if np.sum(mask) >= 2:  # Ensure enough valid data points
                    gene_stabilities_across_datasets[i, j, g], _ = pearsonr(goi_embeddings[mask, i], goi_embeddings[mask, j])

            # Set self-comparison to 1 if gene exists in dataset
            if robust_gene_names[g] in sm_adata_set[i].var_names:
                gene_stabilities_across_datasets[i, i, g] = 1.0

    return gene_stabilities_across_datasets


# Adds a column to the node label dataframe for the minimum gene stability across dataset pairs for each gene
def add_gene_stability_labels(node_label_df, gene_stabilities_across_datasets, robust_gene_names):
    """
    Adds a column to the node label DataFrame for the minimum gene stability across dataset pairs.
    
    Parameters:
    node_label_df (pd.DataFrame): DataFrame containing gene metadata.
    gene_stabilities_across_datasets (np.ndarray): Tensor containing gene stability values.
    robust_gene_names (np.ndarray): Array of gene names corresponding to axis 0 of `gene_embeddings_tensor`.
    
    Returns:
    None: Updates `node_label_df` in-place by adding a 'gene_stabilities' column.
    """

    # Initialize results array with NaNs
    min_gene_stability_arr = np.full(len(node_label_df), np.nan)
    robust_gene_to_index = {gene: idx for idx, gene in enumerate(robust_gene_names)}

    # Iterate through genes in node_label_df and add minimum stability if gene is present in robust_gene_names
    for i, gene in enumerate(node_label_df['name']):
        if gene in robust_gene_names:
            curr_gene_stabilities = gene_stabilities_across_datasets[:,:, robust_gene_to_index[gene]]
            lower_tri_vals = curr_gene_stabilities[np.tril_indices(gene_stabilities_across_datasets.shape[0], k=-1)]
            if np.any(np.isnan(lower_tri_vals)):
                min_gene_stability_arr[i] = -1.0
            else:
                min_gene_stability_arr[i] = np.min(lower_tri_vals)

    # Add new column to node_label_df in place
    node_label_df['gene_stabilities'] = min_gene_stability_arr

