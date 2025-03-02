# Copyright (c) 2025 Chase Holdener
# Licensed under the MIT License. See LICENSE file for details.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import csv
import igraph as ig
import time
import os

from scipy.sparse import lil_matrix, csr_matrix, hstack
from scipy.stats import pearsonr, t


def compute_correlation_matrix(X):
    """
    Computes the pairwise Pearson correlation coefficients between all genes.
    
    Parameters:
    X (numpy.ndarray):
        A (N x G) matrix where rows represent spatial points (N) and columns represent genes (G).
        This matrix should be smoothed gene expression data.
    
    Returns:
    tuple:
        pearsonR_mat (numpy.ndarray): (G x G) matrix where pearsonR_mat[i, j] is the Pearson correlation 
        coefficient between gene i and gene j.
        
        p_val_mat (numpy.ndarray): (G x G) matrix where p_val_mat[i, j] is the p-value associated 
        with the Pearson correlation coefficient (no FDR correction applied).
    """
    print("compute_correlation_matrix running...")
    t0 = time.time()
    
    # get dimensions of X
    N, G = X.shape
    
    # normalize gene surfaces to have mean expression = 0
    means = np.mean(X, axis=0, keepdims=True)
    np.subtract(X, means, out=X)

    # get inner product matrix (G x G)
    I = np.dot(X.T, X)

    # get diagonal of I (1 x G)
    d = np.diag(I)

    # get outer product mat from diagonal (G x G)
    D = np.outer(d, d)

    # get sqrt of D (G x G)
    D_sqrt = np.sqrt(D)

    # calculate Pearson correlation matrix
    pearsonR_mat = I / D_sqrt

    # ensure values of -1 and 1 don't cause division by zero
    eps = 1e-5
    pearsonR_mat = np.clip(pearsonR_mat, -1 + eps, 1 - eps)
    
    # calculate the t-statistics
    t_stat_mat = pearsonR_mat * np.sqrt((N - 2) / (1 - pearsonR_mat**2))
    
    # set diagonal to infinity
    np.fill_diagonal(t_stat_mat, np.inf)
    
    # calculate two-tailed p-values
    p_val_mat = 2 * (1 - t.cdf(np.abs(t_stat_mat), df=N-2))

    # restore original values of X
    np.add(X, means, out=X)

    # runtime report
    print(f"Total runtime for compute_correlation_matrix: {time.time() - t0} seconds.")

    return pearsonR_mat, p_val_mat



def get_correlations_to_GOI(pearsonR_mat, gene_names, GOI, reverse_order=False):
    """
    Retrieves and ranks the correlation of all genes with a specified gene of interest (GOI).
    
    Parameters:
    pearsonR_mat (numpy.ndarray):
        A (G x G) Pearson correlation coefficient matrix.
    gene_names (list, numpy.ndarray, or pd.Index):
        A list or array of gene names corresponding to matrix indices.
    GOI (str):
        The gene of interest for which correlations are ranked.
    reverse_order (bool, optional):
        If True, sorts correlations in ascending order. Defaults to descending order.
    
    Returns:
    numpy.ndarray:
        A (G x 2) array where the first column contains gene names and the second column contains 
        their correlation values with the GOI, sorted based on correlation strength.
    """
    # Validate input data
    assert pearsonR_mat is not None, "pearsonR_mat must be specified."
    assert isinstance(pearsonR_mat, np.ndarray), "pearsonR_mat must be a numpy array."
    assert pearsonR_mat.shape[0] == pearsonR_mat.shape[1], "pearsonR_mat must be a square matrix."
    assert gene_names is not None, "gene_names must be specified."
    assert isinstance(gene_names, (pd.Index, list, tuple, np.ndarray)), "gene_names must be an iterable type."
    assert len(gene_names) == pearsonR_mat.shape[1], "gene_names length must match pearsonR_mat dimensions."
    assert GOI in gene_names, "Error: GOI is not in gene_names."
    
    gene_names = np.array(gene_names)
    goi_num = np.where(gene_names == GOI)[0]
    goi_pearson_R = pearsonR_mat[:, goi_num].flatten()
    
    # Plot histogram of Pearson correlations
    plt.figure(figsize=(6, 4), dpi=300)
    plt.hist(goi_pearson_R, bins=100, color='blue', edgecolor='black')
    plt.xlabel('PCC')
    plt.ylabel('Frequency')
    plt.title(f'Correlation of all genes with {GOI}')
    plt.tight_layout()
    plt.show()
    
    # Sort correlations
    sorted_idx = np.argsort(goi_pearson_R if reverse_order else -goi_pearson_R)
    goi_correlations = np.column_stack((gene_names[sorted_idx], goi_pearson_R[sorted_idx]))
    
    return goi_correlations

