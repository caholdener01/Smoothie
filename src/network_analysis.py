# Copyright (c) 2025 Chase Holdener
# Licensed under the MIT License. See LICENSE file for details.

import pandas as pd
import numpy as np
import random
import csv
import igraph as ig
import time
import os
import matplotlib.pyplot as plt

from scipy.stats import pearsonr, t


def make_spatial_network(pearsonR_mat,
                         gene_names,
                         pcc_cutoff,
                         clustering_power,
                         output_folder,
                         gene_labels_list=None,
                         gene_labels_names=None,
                         save_network=True,
                         save_gene_labels=True,
                         trials=20,
                         random_seed=0):
    """
    Generates a spatial network from a Pearson correlation matrix, performs Infomap clustering, 
    and calculates various network metrics for genes in the dataset.
    
    Parameters:
    -----------
    pearsonR_mat : np.ndarray
        A square matrix of Pearson correlation coefficients between genes.
    
    gene_names : list or pd.Index
        A list of gene names corresponding to the rows/columns of the Pearson correlation matrix.
    
    pcc_cutoff : float
        The threshold for the Pearson correlation coefficient. Only correlations above this value are retained.
    
    clustering_power : float, optional
        A parameter that controls the rescaling of Pearson correlation values. Defaults to 4 if None.

    output_folder : str
        Directory where output files (network and node labels) will be saved.

    gene_labels_list : list of list/tuple/np.ndarray, optional
        A list containing gene set labels. Each item in the list should have the same length as the number of genes.
    
    gene_labels_names : list, optional
        A list of names for the gene sets in `gene_labels_list`.
    
    save_network : bool, optional
        Whether to save the network edge list as a CSV file. Defaults to True.
    
    save_gene_labels : bool, optional
        Whether to save the gene labels as a CSV file. Defaults to True.
    
    trials : int, optional
        Number of trials for the Infomap clustering algorithm. Defaults to 20 if None.
    
    random_seed : int, optional
        Seed for the random number generator used in the Infomap clustering algorithm. Defaults to 0 if None.
    
    Returns:
    --------
    edge_list : list
        List of edges in the format [gene1, gene2, PCC, Rescaled_PCC].
    
    node_label_df : pd.DataFrame
        DataFrame with gene names, community labels, and various network metrics.
    """
    
    ### ASSERTIONS and DEFAULT VARIABLES

    # pearsonR_mat
    assert pearsonR_mat is not None, "pearsonR_mat must be specified."
    assert isinstance(pearsonR_mat, np.ndarray), "pearsonR_mat must be a numpy array."
    assert pearsonR_mat.shape[0] == pearsonR_mat.shape[1], (
        "pearsonR_mat must be square. (pearsonR_mat.shape[0] == pearsonR_mat.shape[1])"
    )
    
    # gene_names
    assert gene_names is not None, "gene_names must be specified. Use sm_adata.var_names."
    assert isinstance(gene_names, (pd.Index, list, set, tuple, np.ndarray)), (
        "gene_names must be a pandas.Index, set, list, tuple, or array."
    )
    assert len(gene_names) == pearsonR_mat.shape[1], (
        "len(gene_names) must equal the column count of pearsonR_mat."
    )
    
    # pcc_cutoff
    assert isinstance(pcc_cutoff, float), "pcc_cutoff must be a floating point number."
    
    # clustering_power
    assert isinstance(clustering_power, (int, float)) and clustering_power >= 1, (
        "clustering_power must be positive. Values between 1 and 6 are reasonable for increasing network modularity."
    )

    # output_folder
    assert isinstance(output_folder, str) and os.path.isdir(output_folder), (
        "Please specify a valid existing directory for output_folder."
    )
    
    # gene_labels_list
    if gene_labels_list is not None:
        assert isinstance(gene_labels_list, list), "gene_labels_list must be a list."
        assert all(isinstance(item, (list, tuple, np.ndarray)) for item in gene_labels_list), (
            "gene_labels_list must be a list containing lists, tuples, or arrays."
        )
        assert all(len(item) == pearsonR_mat.shape[0] for item in gene_labels_list), (
            "Elements of gene_labels_list must have the same length as the number of pearsonR_mat rows."
        )
    
    # gene_labels_names
    if gene_labels_names is not None:
        assert isinstance(gene_labels_names, list), "gene_labels_names must be a list of strings."
        assert len(gene_labels_list) == len(gene_labels_names), (
            "gene_labels_list and gene_labels_names must have the same length."
        )
    
    # trials
    assert isinstance(trials, int) and trials >= 1, "trials must be a positive integer."
    
    
    ### CREATE CORRELATION NETWORK

    gene_names = np.array(gene_names)
    
    # extract lower triangle of pearsonR_mat
    lower_tri_indices = np.tril_indices(pearsonR_mat.shape[0], -1)
    lower_tri_values = pearsonR_mat[lower_tri_indices]
    
    # apply hard thresholding
    valid_indices = np.where(lower_tri_values > pcc_cutoff)[0]
    
    # apply soft thresholding
    rescaled_pearsonR = (lower_tri_values[valid_indices] - pcc_cutoff) / (1.0 - pcc_cutoff)
    rescaled_pearsonR_exp = np.power(rescaled_pearsonR, clustering_power)
    
    # get valid edge indices
    rows, cols = lower_tri_indices[0][valid_indices], lower_tri_indices[1][valid_indices]
    
    # construct edge list: [gene1, gene2, PCC, Rescaled_PCC]
    edge_list = [
        [gene_names[rows[i]], gene_names[cols[i]], pearsonR_mat[rows[i], cols[i]], rescaled_pearsonR_exp[i]]
        for i in range(len(valid_indices))
    ]
    
    # save network edge list
    if save_network:
        output_network_filename = (
            f"{output_folder}/network_PCC{pcc_cutoff}_clustPow{clustering_power}.csv"
        )
        with open(output_network_filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['source', 'target', 'PCC', 'Rescaled_PCC'])  # Write header
            writer.writerows(edge_list)

    
    ### RUN INFOMAP AND CLUSTER GENES INTO MODULES
    
    geneIndexer = pd.Index(gene_names)
    
    # make an array of edges with just the gene indices
    edgeArr_geneNames = np.array(edge_list)[:,0:2]
    edgeArr_geneIndices = np.array(
        (geneIndexer.get_indexer(edgeArr_geneNames[:,0]), geneIndexer.get_indexer(edgeArr_geneNames[:,1]))
    ).T
    
    # construct graph
    g_infomap = ig.Graph(n=len(gene_names), edges=edgeArr_geneIndices)
    g_infomap.es['weight'] = np.array(edge_list)[:, 2].astype(float)  # Original edge weights
    g_infomap.es['rescaled_weight'] = np.array(edge_list)[:, 3].astype(float)  # Monotonically rescaled weights
    
    # run Infomap clustering
    random.seed(random_seed) # set seed for infomap clustering
    infomap = g_infomap.community_infomap(edge_weights=g_infomap.es['rescaled_weight'], trials=trials)
    membership = infomap.membership
    
    # reorder clusters by decreasing cluster size
    community_counts = pd.Series(membership).value_counts()
    sorted_communities = community_counts.sort_values(ascending=False).index
    community_mapping = {old_label: new_label for new_label, old_label in enumerate(sorted_communities, start=1)}
    renumbered_assignments = np.array([community_mapping[label] for label in membership])
    
    # add each metric to the DataFrame
    node_label_df = pd.DataFrame({
        'name': gene_names,
        'community_label': renumbered_assignments,
    })
    
    # calculate some metrics (many more possibilities here!!)
    node_label_df['degree'] = g_infomap.degree()
    node_label_df['weighted_degree'] = g_infomap.strength(weights='weight')
    node_label_df['rescaled_weighted_degree'] = g_infomap.strength(weights='rescaled_weight')
    node_label_df['clustering_coeff'] = g_infomap.transitivity_local_undirected(mode="zero")

    # add other gene set labels to dataframe according to input gene_labels_list
    if gene_labels_list is not None:
        for i, list_i in enumerate(gene_labels_list):
            node_label_df[gene_labels_names[i]] = list_i
            
    # sort by community size and add to node_label_df
    node_label_df.sort_values(by='community_label', inplace=True)
    node_label_df.reset_index(drop=True, inplace=True)

    # save gene labels
    if save_gene_labels:
        output_nodelabels_filename = (
            f"{output_folder}/nodelabels_PCC{pcc_cutoff}_clustPow{clustering_power}.csv"
        )
        node_label_df.to_csv(output_nodelabels_filename, index=False)

    return edge_list, node_label_df


def make_geneset_spatial_network(pearsonR_mat, 
                                 gene_names, 
                                 node_label_df, 
                                 gene_list, 
                                 low_pcc_cutoff, 
                                 output_folder,
                                 intra_geneset_edges_only=False, 
                                 save_network=True, 
                                 save_gene_labels=True):
    """
    Given a clustered network, construct a new network with a lower PCC cutoff for just a provided geneset.

    Parameters:
    -----------
    pearsonR_mat (np.ndarray):  
        Pearson correlation coefficient matrix for genes.  

    gene_names (list):  
        List of gene names corresponding to rows/columns in pearsonR_mat.  

    node_label_df (pd.DataFrame):  
        DataFrame containing gene names and their community labels. 

    gene_list (list):  
        List of genes to retain in the subset network.  

    low_pcc_cutoff (float):  
        Minimum PCC required to include an edge in the new network. 

    output_folder : str
        Directory where output files (network and node labels) will be saved.

    intra_geneset_edges_only (bool):
        A boolean value specifying what edges should be added to the network. Defaults to False.
        - True: Only add edges where both nodes appear in gene_list
        - False: Add edges where at least one node appears in gene_list

    save_network : bool, optional
        Whether to save the network edge list as a CSV file. Defaults to True.
    
    save_gene_labels : bool, optional
        Whether to save the gene labels as a CSV file. Defaults to True.

    Returns:
    --------
    geneset_edge_list (list):  
        List of geneset edges in the format [gene1, gene2, PCC].

    geneset_node_label_df (pd.DataFrame):  
        DataFrame with updated community labels, including weak_community_label.
    """
    gene_names = np.array(gene_names)
    
    # Create a mapping from gene name to index
    gene_index_map = {gene: i for i, gene in enumerate(gene_names)}

    # Filter Pearson matrix for genes in and out of the gene list
    geneset_indices = [gene_index_map[gene] for gene in gene_list if gene in gene_index_map]
    non_geneset_indices = [gene_index_map[gene] for gene in gene_names if gene not in gene_index_map]

    ## Construct edgelist with intra-geneset edges
    sub_mat = pearsonR_mat[np.ix_(geneset_indices, geneset_indices)]

    # extract lower triangle of pearsonR_mat
    lower_tri_indices = np.tril_indices(sub_mat.shape[0], -1)
    lower_tri_values = sub_mat[lower_tri_indices]
    
    # apply hard thresholding
    valid_indices = np.where(lower_tri_values > low_pcc_cutoff)[0]
    
    # get valid edge indices
    rows, cols = lower_tri_indices[0][valid_indices], lower_tri_indices[1][valid_indices]
    
    # construct edge list: [gene1, gene2, PCC]
    geneset_edge_list = [
        [gene_names[geneset_indices[rows[i]]], gene_names[geneset_indices[cols[i]]], sub_mat[rows[i], cols[i]]]
        for i in range(len(valid_indices))
    ]

    ## Optionally, construct edgelist between gene_list members and non-gene_list members
    if not intra_geneset_edges_only:
        
        outer_sub_mat = pearsonR_mat[np.ix_(geneset_indices, non_geneset_indices)]

        # apply hard thresholding
        valid_indices = np.where(outer_sub_mat > low_pcc_cutoff)[0]
        
        # get valid edge indices
        rows, cols = outer_sub_mat[0][valid_indices], outer_sub_mat[1][valid_indices]
        
        # construct edge list: [gene1, gene2, PCC]
        inter_geneset_edge_list = [
            [gene_names[geneset_indices[rows[i]]], gene_names[non_geneset_indices[cols[i]]], outer_sub_mat[rows[i], cols[i]]]
            for i in range(len(valid_indices))
        ]
        # Concatenate edgelists for output
        geneset_edge_list = geneset_edge_list + inter_geneset_edge_list

    ## Construct geneset_node_label_df
    # Filter node_label_df to keep only communities with size >= 2
    filtered_communities = node_label_df.groupby('community_label').filter(lambda x: len(x) >= 2)

    # Identify strongly connected genes
    strong_genes = set(filtered_communities['name'])

    # Initialize weak_community_label column
    filtered_communities['weak_community_label'] = -1

    # Adjust diagonal of PCC matrix to 0.0
    np.fill_diagonal(pearsonR_mat, 0.0)
    
    # Assign weak_community_labels for gene_list genes not in strong_genes
    for gene in gene_list:
        if gene not in strong_genes:
            gene_idx = gene_index_map.get(gene)
            if gene_idx is not None:
                # Get PCC values and find the strongest connected gene in strong_genes
                gene_pcc_values = pearsonR_mat[gene_idx]
                max_index = np.argmax(gene_pcc_values)  # Get index of highest PCC
                highest_pcc_gene = gene_names[max_index]  # Get the corresponding gene name
                highest_pcc_value = gene_pcc_values[max_index]  # Get the highest PCC value
                
                if highest_pcc_gene in strong_genes and highest_pcc_value > low_pcc_cutoff:
                    # Determine weak community label from the strongest connected gene's community label
                    label = filtered_communities.loc[
                        filtered_communities['name'] == highest_pcc_gene, 'community_label'
                    ].values[0]
                    # Add gene_list gene to dataframe
                    new_row = pd.DataFrame([{col: -1 for col in filtered_communities.columns}])
                    new_row.loc[0, ['name', 'weak_community_label']] = [gene, label]
                    filtered_communities = pd.concat([filtered_communities, new_row], ignore_index=True)

    # Adjust diagonal of PCC matrix back to 1.0
    np.fill_diagonal(pearsonR_mat, 1.0)

    # Add column to df for whether gene is a gene_list member
    geneset_member = [gene in gene_list for gene in filtered_communities['name']]
    filtered_communities['geneset_member'] = geneset_member
    # Sort to put gene_list members at the top of output
    geneset_node_label_df = filtered_communities.sort_values(by='geneset_member', ascending=False, kind='stable')

    # save network edge list
    if save_network:
        edgelist_tag = "strict_" if intra_geneset_edges_only else ""
        output_network_filename = (
            f"{output_folder}/{edgelist_tag}geneset_network_low_PCC{low_pcc_cutoff}.csv"
        )
        with open(output_network_filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['source', 'target', 'PCC'])  # Write header
            writer.writerows(geneset_edge_list)
            
    # save gene labels
    if save_gene_labels:
        output_nodelabels_filename = (
            f"{output_folder}/geneset_nodelabels_low_PCC{low_pcc_cutoff}.csv"
        )
        geneset_node_label_df.to_csv(output_nodelabels_filename, index=False)

    return geneset_edge_list, geneset_node_label_df



### Scale-Free Topography Index Functions (adapted from WGCNA - https://pubmed.ncbi.nlm.nih.gov/16646834/)

def compute_sfti(degrees):
    """
    Compute the Scale-Free Topology Index (SFTI) using the R² value 
    of the log-log degree distribution fit to a power-law model.
    """
    if len(degrees) == 0:
        return 0  # No nodes in network

    bins = np.linspace(min(degrees), max(degrees), 100)

    hist, bin_edges = np.histogram(degrees, bins=bins)

    p_k = hist / hist.sum()

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    mask = (bin_centers > 0) & (p_k > 0)

    # Apply log only on valid values
    log_bin_centers = np.log(bin_centers[mask])
    log_p_k = np.log(p_k[mask])
    
    # Compute the correlation
    corr = np.corrcoef(log_bin_centers, log_p_k)[0, 1]

    sfti = corr ** 2

    return sfti

def analyze_sfti_vs_thresholding_params(corr_matrix, cutoffs, powers):
    """
    Iterates over different correlation cutoffs and soft thresholding powers,
    constructs networks, computes SFTI, and returns an array of inputs vs SFTI values.
    """
    sfti_results = []

    for cutoff in cutoffs:
        
        for power in powers:
            
            adj_matrix = (corr_matrix >= cutoff).astype(int)
    
            np.fill_diagonal(adj_matrix, 0)

            # Step 1: Rescale values in the range [cutoff, 1] to [0, 1], avoid division by zero
            corr_matrix_rescaled = np.copy(corr_matrix)
            mask = (corr_matrix_rescaled >= cutoff) & (corr_matrix_rescaled <= 1)
            
            # Apply the rescaling, only where the condition is true
            corr_matrix_rescaled[mask] = (corr_matrix_rescaled[mask] - cutoff) / (1 - cutoff)
        
            # Step 2: Handle NaNs and raise to the power of P (ignoring np.nan values)
            corr_matrix_rescaled = np.where(np.isnan(corr_matrix_rescaled), np.nan, corr_matrix_rescaled ** power)
            
            # Compute weighted degrees
            weighted_degrees = np.sum(corr_matrix_rescaled * adj_matrix, axis=1)
            
            # Filter out weighted degrees of zero or negative
            weighted_degrees_filtered = weighted_degrees[weighted_degrees > 0]
    
            # Compute SFTI
            sfti = compute_sfti(weighted_degrees)
            sfti_results.append((cutoff, power, sfti))

    # make sfti_results an array
    sfti_results = np.array(sfti_results)
    
    # Plot results
    plt.figure(figsize=(6, 4))

    sfti_results = np.array(sfti_results)
    
    for power in powers:
        # Filter data for this power
        sfti_for_power = sfti_results[sfti_results[:, 1] == power]
        
        # Extract cutoff values and SFTI values for this power
        cutoffs_for_power = sfti_for_power[:, 0]
        sfti_values_for_power = sfti_for_power[:, 2]
        
        # Plot the SFTI values for the given power
        plt.plot(cutoffs_for_power, sfti_values_for_power, marker='o', label=f'{power}')
    
    # Add labels and title
    plt.xlabel('Hard Correlation Cutoff')
    plt.ylabel('SFTI (R²)')
    plt.title('SFTI across thresholding inputs')
    plt.legend(title='Soft Power', loc='upper left', bbox_to_anchor=(1, 1))
    plt.grid(True)
    plt.show()

    return sfti_results
