# Copyright (c) 2026 Chase Holdener
# Licensed under the MIT License. See LICENSE file for details.

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import anndata as ad
import os
import warnings
from scipy.stats import pearsonr


def concatenate_smoothed_matrices(sm_adata_list):
    """
    Concatenates smoothed count matrices from multiple datasets, aligning gene columns.
    
    Parameters:
    sm_adata_list (list of AnnData): List of AnnData objects with smoothed count matrices.
    
    Returns:
    tuple: 
        - X_concat (np.ndarray): Concatenated matrix with genes aligned across datasets.
        - gene_names (list of str): List of gene names, each corresponding to a col of X_concat.
    """
    print("concatenate_smoothed_matrices running...")
    
    # Collect unique gene names and sort them
    gene_names = sorted(set.union(*(set(sm_adata.var_names) for sm_adata in sm_adata_list)))
    gene_names_to_index = {gene: idx for idx, gene in enumerate(gene_names)}

    # Determine total obs and precompute obs index ranges
    obs_counts = np.array([sm_adata.shape[0] for sm_adata in sm_adata_list])
    sm_adata_obs_ranges = np.concatenate(([0], np.cumsum(obs_counts)))

    # Preallocate matrix
    X_concat = np.zeros((sm_adata_obs_ranges[-1], len(gene_names)), dtype=np.float32)

    # Fill the matrix
    for i, sm_adata in enumerate(sm_adata_list):
        row_start, row_end = sm_adata_obs_ranges[i], sm_adata_obs_ranges[i+1]
        col_indices = [gene_names_to_index[gene] for gene in sm_adata.var_names]
        X_concat[row_start:row_end, col_indices] = sm_adata.X

    return X_concat, gene_names


def align_gene_correlation_matrices(sm_adata_list, pearsonR_mat_list):
    """
    Aligns gene correlation matrices such that each row and column index corresponds to the same gene across all matrices.
    
    Parameters:
    sm_adata_list (list of AnnData): List of AnnData objects with smoothed count matrices.
    pearsonR_mat_list (list of np.ndarray): List of Pearson correlation matrices corresponding to each dataset.
    
    Returns:
    aligned_pearsonR_tensor (np.ndarray): A 3D tensor where each slice along the third axis is an aligned Pearson correlation matrix.
    """
    
    # Collect unique gene names and sort them
    gene_names = sorted(set.union(*(set(sm_adata.var_names) for sm_adata in sm_adata_list)))
    gene_names_to_index = {gene: idx for idx, gene in enumerate(gene_names)}
    
    # Initialize a 3D tensor of pearsonR matrices for each dataset as np.nan
    aligned_pearsonR_tensor = np.full((len(gene_names), len(gene_names), len(sm_adata_list)), 
                                      np.nan, dtype=np.float32)
    
    # Fill in the tensor for each dataset
    for d, sm_adata in enumerate(sm_adata_list):
        # Get indices of genes that exist in gene_names
        mask = np.array([gene in gene_names_to_index for gene in sm_adata.var_names])
        valid_genes = np.array(sm_adata.var_names)[mask]
        
        # Get corresponding indices
        idx_in_global = np.array([gene_names_to_index[gene] for gene in valid_genes])
        idx_in_local = np.where(mask)[0]

        # Extract the submatrix from the current dataset
        submatrix = pearsonR_mat_list[d][idx_in_local[:, None], idx_in_local]
        
        # Assign values using NumPy advanced indexing
        aligned_pearsonR_tensor[idx_in_global[:, None], idx_in_global, d] = submatrix 

    return aligned_pearsonR_tensor


def create_gene_embeddings(aligned_pearsonR_tensor, sm_adata_list, pcc_cutoff, node_label_df=None, E_max=25, seed=0):
    """
    Extracts gene embeddings by filtering aligned Pearson correlation matrices for spatially informative genes  
    and optionally balancing embedding features based on gene module assignments.
    
    Parameters:
    -----------
    aligned_pearsonR_tensor (np.ndarray):  
        A 3D tensor containing aligned Pearson correlation matrices across datasets.  
    
    sm_adata_list (list of AnnData):  
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
            - Axis 2: Spatial datasets, in the same order as `sm_adata_list`.  
        
        robust_gene_names (np.ndarray):  
            Array of gene names corresponding to the first axis (embedding vectors) of `gene_embeddings_tensor`.  
    """
    
    # Collect unique gene names and sort them
    gene_names = sorted(set.union(*(set(sm_adata.var_names) for sm_adata in sm_adata_list)))
    
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
    
    # Set correlations >= 0.999 or <= -0.999 to NaN before Fisher Z-transform
    aligned_pearsonR_tensor_robust[np.abs(aligned_pearsonR_tensor_robust) >= 0.999] = np.nan
    
    # Convert to Fisher's Z-Scores
    Zscore_aligned_pearsonR_tensor_robust = 0.5 * np.log((1 + aligned_pearsonR_tensor_robust) / (1 - aligned_pearsonR_tensor_robust))

    # Downsample large gene module embedding features if modules are provided
    if node_label_df is not None and E_max is not None:
        def sample_indices(group):
            return np.random.choice(group.index, min(len(group), E_max), replace=False)

        np.random.seed(seed)

        try:
            # Attempt the modern Pandas 2.2.0+ syntax
            sampled_indices = node_label_df.groupby('module_label', group_keys=False).apply(
                sample_indices, include_groups=False
            ).explode().astype(int)
    
        except TypeError:
            # Fall back to Pandas < 2.2.0 syntax if include_groups is rejected
            sampled_indices = node_label_df.groupby('module_label', group_keys=False).apply(
                sample_indices
            ).explode().astype(int)
        
        # Create a mask for selected genes
        selected_genes = set(node_label_df.loc[sampled_indices, 'name'])
        downsampling_mask = np.array([gene in selected_genes for gene in robust_gene_names])

        # Create gene embeddings with downsampling
        gene_embeddings_tensor = Zscore_aligned_pearsonR_tensor_robust[:, downsampling_mask, :]

    else: # Take gene embeddings without module-based feature downsampling 
        gene_embeddings_tensor = Zscore_aligned_pearsonR_tensor_robust

    # Return gene embeddings and the gene names for each row vector
    return gene_embeddings_tensor, robust_gene_names


def get_gene_stabilities_across_datasets(gene_embeddings_tensor, robust_gene_names, sm_adata_list):
    """
    Computes stability of same gene embeddings across dataset pairs using Pearson correlation.
    
    Parameters:
    -----------
    gene_embeddings_tensor (np.ndarray): 
        Tensor of gene embeddings output from function create_gene_embeddings.
    
    robust_gene_names (np.ndarray):
        Array of gene names corresponding to axis 0 of `gene_embeddings_tensor`.
    
    sm_adata_list (list of AnnData):
        List of AnnData objects containing smoothed count matrices.
    
    Returns:
    --------
    gene_stabilities_across_datasets (np.ndarray):
        A tensor containing stability values for each gene across dataset pairs.
        - Axis 0: Spatial datasets, in the same order as `sm_adata_list`
        - Axis 1: Spatial datasets, in the same order as `sm_adata_list`
        - Axis 2: Genes, in same order as robust_gene_names.

        gene_stabilities_across_datasets[i,j,k] is the stability of gene
        robust_gene_names[k] between dataset sm_adata_list[i] and dataset sm_adata_list[j] !
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
            if robust_gene_names[g] in sm_adata_list[i].var_names:
                gene_stabilities_across_datasets[i, i, g] = 1.0

    return gene_stabilities_across_datasets


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


def run_second_order_correlation_analysis(sm_adata_list, pcc_cutoff, node_label_df=None, E_max=25, 
                                         output_folder=None, seed=0):
    """
    Comprehensive function for second-order correlation gene embedding analysis across multiple samples.
    Workflow: 1) aligns correlation matrices, 2) creates gene embeddings, 
    3) computes stability metrics, and 4) adds stability values to node_label_df.
    
    Parameters:
    -----------
    sm_adata_list : list of AnnData
        List of AnnData objects containing smoothed count matrices from different spatial transcriptomics samples.
        Each AnnData should have gene expression data in .X and gene names in .var_names.
    
    pcc_cutoff : float
        Minimum Pearson correlation coefficient threshold. Genes must exceed this correlation value in at least
        one dataset to be included in the embedding analysis. Typical values range from 0.4 to 0.5.
    
    node_label_df : pd.DataFrame, optional
        DataFrame containing multi-dataset gene module assignments with columns 'name' (gene names) and 
        'module_label' (module IDs). If provided:
        - Embedding features will be balanced across modules
        - Stability metrics will be automatically added as 'gene_stabilities' column
        If None, all genes above pcc_cutoff are used as features without module-based downsampling. 
        Default is None.

    E_max : int
        Maximum number of genes per module in node_label_df to use for gene embedding features. 
        This parameter controls the dimensionality of gene embedding vectors and helps balance 
        representation across gene modules. Typical values: 15-50.
    
    output_folder : str, optional
        Path to folder where results will be saved. 
        If None, files are not saved. Default is None.
    
    seed : int, optional
        Random seed for reproducibility when downsampling embedding features. Default is 0.
    
    Returns:
    --------
    dict
        A dictionary containing all analysis results:
        
        - 'gene_embeddings_tensor' (np.ndarray): 3D tensor of gene embeddings with shape 
          (n_genes, n_features, n_datasets) where each gene's embedding is a vector of Fisher Z-transformed
          correlation values.
        
        - 'robust_gene_names' (np.ndarray): Array of gene names corresponding to rows in gene_embeddings_tensor.
          These are the genes that passed the pcc_cutoff threshold in at least one dataset in sm_adata_list.
        
        - 'gene_stabilities' (np.ndarray): 3D tensor with shape (n_datasets, n_datasets, n_genes) containing
          pairwise correlation values between gene embeddings across datasets. Values range from -1 to 1,
          where higher values indicate more stable/conserved co-expression patterns.
        
        - 'aligned_pearsonR_tensor' (np.ndarray): 3D tensor with shape (n_all_genes, n_all_genes, n_datasets)
          containing aligned Pearson correlation matrices across all datasets.
        
        - 'node_label_df' (pd.DataFrame or None): If node_label_df was provided as input, returns it with
          added 'gene_stabilities' column. Otherwise returns None.
    
    Notes:
    ------
    - Gene embeddings are Fisher Z-transformed correlation profiles, making them suitable for Pearson correlation
      comparison across datasets.
    - The stability metric quantifies how consistently genes co-vary with the same partners across different
      biological samples or experimental conditions.
    - Genes not present in robust_gene_names will have NaN stability values.
    - Genes present in robust_gene_names but with missing data in at least one dataset in sm_adata_list 
      will have stability = -1.0.
    """
    
    from .spatial_correlation import compute_correlation_matrix
    import os
    
    print("Starting second-order correlation analysis pipeline...")
    print(f"Number of datasets: {len(sm_adata_list)}")
    
    # Step 1: Compute correlation matrices for each dataset
    print("\n" + "=" * 60 + "\n[Step 1/5] Computing Pearson correlation matrices for each dataset...")
    pearsonR_mat_list = []
    for i, sm_adata in enumerate(sm_adata_list):
        print(f"\nComputing correlations for dataset {i+1}/{len(sm_adata_list)} ({sm_adata.n_vars} genes, {sm_adata.n_obs} spots)...")
        pearsonR_mat, _ = compute_correlation_matrix(sm_adata.X)
        pearsonR_mat_list.append(pearsonR_mat)
    
    # Step 2: Align correlation matrices across datasets
    print("\n" + "=" * 60 + "\n[Step 2/5] Aligning gene correlation matrices across datasets...")
    aligned_pearsonR_tensor = align_gene_correlation_matrices(sm_adata_list, pearsonR_mat_list)
    print(f"Aligned tensor shape: {aligned_pearsonR_tensor.shape}")
    
    # Step 3: Create gene embeddings
    print("\n" + "=" * 60 + f"\n[Step 3/5] Creating gene embeddings (PCC cutoff={pcc_cutoff}, E_max={E_max})...")
    
    gene_embeddings_tensor, robust_gene_names = create_gene_embeddings(
        aligned_pearsonR_tensor=aligned_pearsonR_tensor,
        sm_adata_list=sm_adata_list,
        pcc_cutoff=pcc_cutoff,
        node_label_df=node_label_df,
        E_max=E_max,
        seed=seed
    )
    print(f"Gene embeddings shape: {gene_embeddings_tensor.shape}")
    
    # Step 4: Compute gene stabilities across datasets
    print("\n" + "=" * 60 + "\n[Step 4/5] Computing gene stability metrics across dataset pairs...")
    gene_stabilities = get_gene_stabilities_across_datasets(
        gene_embeddings_tensor=gene_embeddings_tensor,
        robust_gene_names=robust_gene_names,
        sm_adata_list=sm_adata_list
    )
    print(f"Stability tensor shape: {gene_stabilities.shape}")

    # Add stability labels to node_label_df if provided
    if node_label_df is not None:
        add_gene_stability_labels(node_label_df, gene_stabilities, robust_gene_names)
    
    # Step 5: Save results if output folder specified
    if output_folder is not None:
        print("\n" + "=" * 60 + f"\n[Step 5/5] Saving results to {output_folder}...")
        os.makedirs(output_folder, exist_ok=True)
        
        # Save tensors
        np.save(os.path.join(output_folder, 'gene_embeddings_tensor.npy'), gene_embeddings_tensor)
        np.save(os.path.join(output_folder, 'robust_gene_names.npy'), robust_gene_names)
        np.save(os.path.join(output_folder, 'gene_stabilities_tensor.npy'), gene_stabilities)
        
        # Save node_label_df with stabilities if provided
        if node_label_df is not None:
            output_path = os.path.join(output_folder, 'nodelabels_with_stabilities.csv')
            node_label_df.to_csv(output_path, index=False)
    else:
        print("\n" + "=" * 60 + "\n[Step 5/5] Skipping file output (no output_folder specified)")
    
    # Compile results
    results = {
        'gene_embeddings_tensor': gene_embeddings_tensor,
        'robust_gene_names': robust_gene_names,
        'gene_stabilities': gene_stabilities,
        'aligned_pearsonR_tensor': aligned_pearsonR_tensor,
        'node_label_df': node_label_df if node_label_df is not None else None
    }
    
    print("Second-order correlation analysis complete.")
    
    return results



def plot_dataset_pair_stabilities(gene_stabilities, dataset_names, output_folder=None, 
                                  figsize=None, fontsize=5, dpi=300, file_format='png'):
    """
    Visualize distribution of gene stability between all dataset pairs using violin plots.
    
    Parameters
    ----------
    gene_stabilities : np.ndarray
        Gene stabilities tensor of shape (n_datasets, n_datasets, n_genes)
        from run_second_order_correlation_analysis().
    dataset_names : list of str
        Names of the datasets for axis labels.
    output_folder : str, optional
        Directory to save the plot. If None, does not save.
    figsize : tuple, optional
        Figure size in inches (width, height). If None, auto-calculated based on number of pairs.
    fontsize : int, optional
        Font size for labels and titles. Default is 5.
    dpi : int, optional
        Resolution for saved figure. Default is 300.
    file_format : str, optional
        Format to save the figure: 'png' or 'pdf'. Default is 'png'.
    """
    
    n_datasets = gene_stabilities.shape[0]
    
    # Initialize a list to store the dataframes
    rows_list = []
    
    # Collect data for each dataset pair (lower triangle only - unique pairs)
    for i in range(n_datasets):
        for j in range(i):  # j < i to get lower triangle
            pair_name = f'{dataset_names[j]}+{dataset_names[i]}'
            stabilities = gene_stabilities[i, j, :]
            
            # Filter out NaN values
            valid_stabilities = stabilities[~np.isnan(stabilities)]
            
            # Only append if there is valid data to avoid empty DF issues
            if len(valid_stabilities) > 0:
                new_rows = pd.DataFrame({
                    'dataset_comparison': [pair_name] * len(valid_stabilities),
                    'gene_stability': valid_stabilities
                })
                rows_list.append(new_rows)
    
    # Concatenate all rows once at the end
    if rows_list:
        stability_df = pd.concat(rows_list, ignore_index=True)
    else:
        # Fallback for empty data
        stability_df = pd.DataFrame(columns=['dataset_comparison', 'gene_stability'])

    # Calculate figure size based on number of pairs
    n_pairs = n_datasets * (n_datasets - 1) // 2
    if figsize is None:
        # Avoid division by zero if n_pairs is 0
        width_factor = float(n_pairs) / 3.2 if n_pairs > 0 else 1.0
        figsize = (width_factor, 1)
    
    # Create violin plot
    plt.figure(figsize=figsize)
    
    # Check if we have data to plot
    if not stability_df.empty:
        sns.violinplot(
            x='dataset_comparison',
            y='gene_stability',
            data=stability_df,
            linewidth=0.6,
            width=0.6
        )
        
        # Customize violin appearance
        for violin in plt.gca().collections:
            violin.set_facecolor('lightgray')
            violin.set_edgecolor('black')
    
    # Formatting
    plt.title('Gene Stability Distributions', fontsize=fontsize, pad=4)
    plt.xlabel('')
    plt.ylabel('Pattern Stability', fontsize=fontsize)
    plt.grid(False)
    plt.xticks(rotation=30, fontsize=fontsize, ha='right')
    plt.yticks([-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1], fontsize=fontsize)
    plt.gca().tick_params(axis='x', labelsize=fontsize, pad=0)
    plt.gca().tick_params(axis='y', labelsize=fontsize, pad=1)
    plt.gca().yaxis.labelpad = 0
    
    # Save plot
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)
        if file_format not in ['png', 'pdf']:
            print(f"Warning: file_format '{file_format}' not recognized. Defaulting to 'png'.")
            file_format = 'png'
            
        filename = os.path.join(output_folder, f'pattern_stability_distributions_pairs.{file_format}')
        plt.savefig(filename, format=file_format, bbox_inches='tight', dpi=dpi, pad_inches=0)
    
    plt.show()



def plot_gene_stability_distribution(gene_stabilities, robust_gene_names, node_label_df=None, output_folder=None,
                                     figsize=(3,2), fontsize=6, bins=50, dpi=300, file_format='png'):
    """
    Plot distribution of minimum gene stability across all dataset pairs.
    
    Parameters
    ----------
    gene_stabilities : np.ndarray
        Gene stabilities tensor of shape (n_datasets, n_datasets, n_genes)
        from run_second_order_correlation_analysis().
    robust_gene_names : np.ndarray
        Array of gene names corresponding to the third dimension of gene_stabilities.
    node_label_df : pd.DataFrame, optional
        Node labels DataFrame with 'name' and 'module_label' columns. If provided,
        module labels will be included in output.
    output_folder : str, optional
        Directory to save the plot. If None, does not save.
    figsize : tuple, optional
        Figure size in inches (width, height). Default is (2, 2).
    fontsize : int, optional
        Font size for axes labels. Default is 6.
    bins : int, optional
        Number of histogram bins. Default is 100.
    dpi : int, optional
        Resolution for saved figure. Default is 300.
    file_format : str, optional
        Format to save the figure: 'png' or 'pdf'. Default is 'png'.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['gene', 'module', 'min_stability'] sorted from most to least stable.
        If node_label_df not provided, 'module' column will be NaN.
    """
    
    n_datasets = gene_stabilities.shape[0]
    n_genes = gene_stabilities.shape[2]
    
    # Initialize with -1.0 as default (indicates missing data)
    min_stabilities = np.full(n_genes, -1.0)
    
    # For each gene, compute minimum stability across all dataset pairs
    for g in range(n_genes):
        # Extract lower triangle (unique pairwise comparisons)
        curr_gene_stabilities = gene_stabilities[:, :, g][
            np.tril_indices(n_datasets, k=-1)
        ]
        
        # Skip if any NaN values (missing data in some comparison)
        if np.any(np.isnan(curr_gene_stabilities)):
            continue
        else:
            min_stabilities[g] = np.min(curr_gene_stabilities)

    with plt.rc_context({
        'font.size': fontsize,
        'axes.titlesize': fontsize,
        'axes.labelsize': fontsize,
        'xtick.labelsize': fontsize,
        'ytick.labelsize': fontsize
    }):
        # Plot histogram (including -1 values)
        plt.figure(figsize=figsize)
        plt.hist(min_stabilities, bins=bins, color='blue', edgecolor='black')
        plt.xlabel('Gene pattern stability across all datasets', fontsize=fontsize)
        plt.ylabel('Frequency', fontsize=fontsize)
        plt.gca().tick_params(axis='x', labelsize=fontsize, pad=0)
        plt.gca().tick_params(axis='y', labelsize=fontsize, pad=1)
        plt.title('Distribution of Minimum Gene Stability', fontsize=fontsize)
        plt.grid(True, alpha=0.3, axis='y')
        plt.tight_layout()
        
        if output_folder:
            os.makedirs(output_folder, exist_ok=True)
            if file_format not in ['png', 'pdf']:
                print(f"Warning: file_format '{file_format}' not recognized. Defaulting to 'png'.")
                file_format = 'png'
                
            filename = os.path.join(output_folder, f'gene_stability_distribution.{file_format}')
            plt.savefig(filename, format=file_format, bbox_inches='tight', dpi=dpi)
    
        plt.show()
    
    # Create DataFrame sorted from most to least stable
    stability_df = pd.DataFrame({
        'gene': robust_gene_names,
        'min_stability': min_stabilities
    })
    
    # Add module labels if node_label_df provided
    if node_label_df is not None:
        # Create mapping from gene name to module
        gene_to_module = dict(zip(node_label_df['name'], node_label_df['module_label']))
        stability_df['module'] = stability_df['gene'].map(gene_to_module)
    else:
        stability_df['module'] = np.nan
    
    # Reorder columns and sort
    stability_df = stability_df[['gene', 'module', 'min_stability']]
    stability_df = stability_df.sort_values('min_stability', ascending=False).reset_index(drop=True)
    
    return stability_df



def plot_gene_stability(dataset_names, gene_name, gene_stabilities, robust_gene_names, output_folder=None,
    figsize=(1.5, 1.5), fontsize=5, dpi=300, file_format='png', x_ticks=True, cbar=True):
    """
    Plot compact stability heatmap for a specific gene showing pairwise dataset correlations.

    Parameters
    ----------
    dataset_names : list of str
        Names of the datasets for axis labels.
    gene_name : str
        Name of the gene to visualize.
    gene_stabilities : np.ndarray
        Gene stabilities tensor of shape (n_datasets, n_datasets, n_genes).
    robust_gene_names : np.ndarray
        Array of gene names corresponding to the third dimension of gene_stabilities.
    output_folder : str, optional
        Directory to save the plot. If None, does not save.
    figsize : tuple, optional
        Figure size in inches (width, height). Default is (1, 1).
    fontsize : int, optional
        Font size for title and tick labels. Default is 5.
    dpi : int, optional
        Resolution for saved figure. Default is 300.
    file_format : str, optional
        Format to save the figure: 'png' or 'pdf'. Default is 'png'.
    x_ticks : bool, optional
        Whether to display x-axis dataset labels. Default is True.
    cbar : bool, optional
        Whether to display a colorbar with label. Default is True.
    """

    # Check if gene exists
    if gene_name not in robust_gene_names:
        print(f"Error: Gene '{gene_name}' not found in robust gene set")
        print(f"Available genes: {len(robust_gene_names)}")
        return

    # Get gene index
    gene_idx = np.where(robust_gene_names == gene_name)[0][0]

    # Extract stability matrix
    gene_stab_matrix = gene_stabilities[:, :, gene_idx]

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(
        gene_stab_matrix,
        cmap='viridis',
        interpolation='nearest',
        vmin=-0.5,
        vmax=1.0
    )

    ax.set_title(gene_name, fontsize=fontsize, pad=2)

    if x_ticks:
        ax.set_xticks(range(len(dataset_names)))
        ax.set_xticklabels(dataset_names, fontsize=fontsize, rotation=45, ha='right')
    else:
        ax.set_xticks([])

    # Y ticks
    ax.set_yticks(range(len(dataset_names)))
    ax.set_yticklabels(dataset_names, fontsize=fontsize)

    ax.tick_params(axis='both', labelsize=fontsize, pad=1)
    ax.grid(False)

    if cbar:
        cbar_obj = fig.colorbar(im, ax=ax, fraction=0.07, pad=0.05, shrink=0.7)
        cbar_obj.set_label("Gene stability", fontsize=fontsize)
        cbar_obj.ax.tick_params(labelsize=fontsize)

    plt.tight_layout(pad=0.2)

    # Save if requested
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

        if file_format not in ['png', 'pdf']:
            print(f"Warning: file_format '{file_format}' not recognized. Defaulting to 'png'.")
            file_format = 'png'

        safe_gene_name = gene_name.replace('/', '_').replace('\\', '_')
        filename = os.path.join(
            output_folder,
            f'{safe_gene_name}_stability_plot.{file_format}'
        )

        plt.savefig(filename, format=file_format, bbox_inches='tight', dpi=dpi)

    plt.show()



def plot_module_stability(module_label, gene_stabilities, robust_gene_names, node_label_df, dataset_names,
    output_folder=None, figsize=(1.5, 1.5), fontsize=5, dpi=300, file_format='png', x_ticks=False, cbar=True
):
    """
    Plot compact stability heatmap for a module showing average pairwise dataset correlations.

    Parameters
    ----------
    dataset_names : list of str
        Names of the datasets for axis labels.
    gene_name : str
        Name of the gene to visualize.
    gene_stabilities : np.ndarray
        Gene stabilities tensor of shape (n_datasets, n_datasets, n_genes).
    robust_gene_names : np.ndarray
        Array of gene names corresponding to the third dimension of gene_stabilities.
    output_folder : str, optional
        Directory to save the plot. If None, does not save.
    figsize : tuple, optional
        Figure size in inches (width, height). Default is (0.5, 0.5).
    fontsize : int, optional
        Font size for title and labels. Default is 5.
    dpi : int, optional
        Resolution for saved figure. Default is 300.
    file_format : str, optional
        Format to save the figure: 'png' or 'pdf'. Default is 'png'.
    x_ticks : bool, optional
        Whether to display x-axis dataset labels. Default is True.
    cbar : bool, optional
        Whether to display a colorbar with label. Default is True.
    """

    # Get genes in this module
    module_genes = node_label_df[
        node_label_df['module_label'] == module_label]['name'].values

    if len(module_genes) == 0:
        print(f"Error: Module {module_label} not found or has no genes")
        return

    # Find indices in robust gene set
    mask = np.isin(robust_gene_names, module_genes)
    gene_indices = np.where(mask)[0]

    if len(gene_indices) == 0:
        print(f"Error: No genes from module {module_label} found in robust gene set")
        return

    print(f"Module {module_label}: {len(gene_indices)} genes")

    # Extract stability matrices
    module_stab_matrices = gene_stabilities[:, :, gene_indices]

    # Compute mean across genes
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        avg_module_stab_matrix = np.nanmean(module_stab_matrices, axis=2)

    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad(color='none')

    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(avg_module_stab_matrix, cmap=cmap, interpolation='nearest',
                   vmin=-0.5, vmax=1.0)

    ax.set_title(f'Module {module_label}', fontsize=fontsize, pad=2)
    if x_ticks:
        ax.set_xticks(range(len(dataset_names)))
        ax.set_xticklabels(dataset_names, fontsize=fontsize, rotation=90)
    else:
        ax.set_xticks([])
    ax.set_yticks(range(len(dataset_names)))
    ax.set_yticklabels(dataset_names, fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize, pad=1, length=0)
    ax.grid(False)

    if cbar:
        cbar_obj = fig.colorbar(im, ax=ax, fraction=0.07, pad=0.05, shrink=0.7)
        cbar_obj.set_label("Module stability", fontsize=fontsize)
        cbar_obj.ax.tick_params(labelsize=fontsize)

    plt.tight_layout(pad=0.2)

    # Save
    if output_folder:
        os.makedirs(output_folder, exist_ok=True)

        if file_format not in ['png', 'pdf']:
            print(f"Warning: file_format '{file_format}' not recognized. Defaulting to 'png'.")
            file_format = 'png'

        filename = os.path.join(
            output_folder,
            f'module_{module_label}_stability_plot.{file_format}'
        )

        plt.savefig(filename, format=file_format, bbox_inches='tight', dpi=dpi,
                    pad_inches=0, transparent=True)

    plt.show()