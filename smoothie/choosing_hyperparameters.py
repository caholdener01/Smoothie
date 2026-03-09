# Copyright (c) 2026 Chase Holdener
# Licensed under the MIT License. See LICENSE file for details.

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import sys

from scipy import ndimage

from .gaussian_smoothing import run_gridbased_smoothing
from .network_analysis import make_spatial_network


def select_clustering_params(
    gene_names,
    pearsonR_mat,
    permuted_pcc_999=None,
    output_folder=None,
    pcc_cutoffs=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
    clustering_powers=[1, 3, 5, 7, 9],
    min_genes_for_module=3,
    infomap_clustering_trials=1,
    full_metrics=False
):

    """
    Evaluates combinations of PCC cutoff and clustering power hyperparameters for spatial 
    network construction, generating diagnostic plots of clustering quality metrics to 
    guide hyperparameter selection.

    Parameters:
    -----------
    gene_names : list, np.ndarray, or pd.Index
        Gene names corresponding to the rows/columns of pearsonR_mat.

    pearsonR_mat : np.ndarray
        A square matrix of Pearson correlation coefficients between genes.

    permuted_pcc_999 : float, optional
        The 99.9th percentile PCC from a shuffled-data null distribution 
        (output of compute_shuffled_correlation_percentiles). If provided, 
        this value is appended to pcc_cutoffs as an additional cutoff to evaluate.
        Default is None.

    output_folder : str, optional
        Directory where output plots and CSV results will be saved.
        If None, plots are displayed but not saved. Default is None.

    pcc_cutoffs : list of float, optional
        List of PCC hard-threshold cutoff values to evaluate.
        Network clustering takes significantly longer for lower pcc_cutoff values.
        Default is [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9].

    clustering_powers : list of int or float, optional
        List of soft-power exponents to evaluate for rescaling PCC edge weights.
        Default is [1, 3, 5, 7, 9].

    min_genes_for_module : int, optional
        Minimum number of genes required for a cluster to be counted as a module.
        Default is 3.

    infomap_clustering_trials : int, optional
        Number of Infomap algorithm trials per hyperparameter combination.
        More trials improve stability of cluster assignments. Default is 1 for runtime efficiency.

    full_metrics : bool, optional
        If True, plots all available metrics (mean_gene_margin, fraction_margin_positive, 
        modularity, n_clusters, n_genes_included). If False, plots only 
        (mean_gene_margin, n_clusters, n_genes_included). Default is False.

    Returns:
    --------
    None
        Displays diagnostic plots and optionally saves them along with a 
        results CSV to output_folder.
    """

    ## Two helper functions--------------------------------------------------------------------

    def get_clustering_stats(adj, node_label_modules_df, pcc_cutoff, clustering_power, min_genes_for_module):
        """
        Compute per-node margin scores and cluster-level cohesion.
    
        Node margin: Fisher z-transformed intra-cluster vs best competing cluster.
        Cluster cohesion: average intra-cluster weight scaled by cluster size.
    
        Parameters
        ----------
        adj : np.ndarray or scipy.sparse matrix
            Square adjacency matrix (n x n), filtered to same genes as node_label_modules_df.
        node_label_modules_df : pd.DataFrame
            DataFrame with columns:
                - 'name': node (gene) names in same order as adj
                - 'module_label': cluster assignment
    
        Returns
        -------
        margin_scores : dict
            Node-wise margin scores {node_name: score}
            Akin to Silhouette scores: more positive score means high confidence in cluster assignment,
                                       0 means node is between two clusters,
                                       more negative score means node may belong to a different cluster.
        positive_margin : list of booleans
            Denotes whether the margin score is positive for each node.
        frac_positive : float
            Fraction of nodes with margin > 0.
        modularity : float
            The weighted modularity of the network.
        n_clusters : int
            Number of unique clusters.
        n_nodes : int
            Number of nodes in the analysis.
        """
        
        adj = (adj - pcc_cutoff) / (1 - pcc_cutoff)
        adj = np.clip(adj, a_min=0, a_max=1)
        adj = np.power(adj, clustering_power) 
    
        node_names = node_label_modules_df["name"].tolist()
        cluster_map = dict(zip(node_names, node_label_modules_df["module_label"]))
        name_to_idx = {name: idx for idx, name in enumerate(node_names)}
        
        margin_scores = {}
        clusters = {}
        # collect genes per cluster
        for name, cl in cluster_map.items():
            clusters.setdefault(cl, []).append(name)
        
        # per-node margins
        for i, name in enumerate(node_names):
            ci = cluster_map[name]
            row = adj[i, :]
        
            cluster_sums = {}
            cluster_counts = {}
            
            for j, w in enumerate(row):
                if i == j:
                    continue
                if w > 0:  # only count positive weights
                    cj = cluster_map[node_names[j]]
                    cluster_sums[cj] = cluster_sums.get(cj, 0.0) + w
                    cluster_counts[cj] = cluster_counts.get(cj, 0) + 1
        
            # average intra-cluster
            A = cluster_sums.get(ci, 0.0) / float(cluster_counts.get(ci, 1))
        
            # best competing cluster
            best_avg = 0.0
            for c, s in cluster_sums.items():
                if c == ci:
                    continue
                count = cluster_counts.get(c, 0)
                if count > 0:
                    avg = s / float(count)
                    if avg > best_avg:
                        best_avg = avg
            B = best_avg
        
            # this is positive when avg intra cluster is larger than avg inter cluster w/ 2nd most similar cluster
            margin_scores[name] = (A - B) / (np.max((np.max((A, B)), 1e-12)))
        
        positive_margin = [m > 0 for m in margin_scores.values()]
        frac_positive = np.mean([m > 0 for m in margin_scores.values()])
    
        # compute weighted modularity
        np.fill_diagonal(adj, 0.0)
        
        m = adj.sum() / 2.0
        k = adj.sum(axis=1)
        modularity = 0.0
        for cl, genes in clusters.items():
            indices = [name_to_idx[g] for g in genes]
            sub_adj = adj[np.ix_(indices, indices)]
            W_c = sub_adj.sum()  # total weight of edges inside cluster
            K_c = k[indices].sum()
            modularity += W_c - (K_c**2) / (2*m)
        modularity /= (2*m)
        
        np.fill_diagonal(adj, 0.99999)

        # convert margin scores to arr
        margin_scores = np.array(list(margin_scores.values()))
        
        # extra metrics
        n_clusters = sum(len(genes) >= min_genes_for_module for genes in clusters.values())
        n_nodes = len(node_names)
        
        return (
            margin_scores,
            positive_margin,
            frac_positive,
            modularity,
            n_nodes,
            n_clusters
        )
    
        
    def plot_cutoff_grid(
        df,
        cutoff_col="pcc_cutoff",
        power_col="clustering_power",
        output_folder=None,
        filename="clustering_plots.png"
    ):
        """
        Make a grid of plots: PCC cutoff (x) vs metrics (y), with clustering_power as multiple lines.
        Saves to output_folder if specified.
    
        Parameters
        ----------
        df : pd.DataFrame
            Must contain:
                - cutoff_col (x-axis)
                - power_col (hue, different lines)
                - other numeric columns = metrics to plot (each one gets its own subplot)
        cutoff_col : str
            Column name for PCC cutoff (x-axis).
        power_col : str
            Column name for clustering power (line grouping).
        output_folder : str or None
            Directory to save the figure. If None, the plot is just shown.
        filename : str
            File name for the saved figure (default "cutoff_grid.png").
        """
        # Identify metric columns (everything except cutoff and power)
        if full_metrics:
            metrics = ["mean_gene_margin", "fraction_margin_positive", "modularity", "n_clusters", "n_genes_included"]
        else:
            metrics = ["mean_gene_margin", "n_clusters", "n_genes_included"]
    
        n_metrics = len(metrics)
        ncols = 3
        nrows = (n_metrics + ncols - 1) // ncols
    
        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows), squeeze=False)
    
        for ax, metric in zip(axes.flat, metrics):

            sns.lineplot(
                data=df,
                x=cutoff_col,
                y=metric,
                hue=power_col,
                marker="o",
                ax=ax
            )
            ax.set_title(metric)
            if metric == "n_clusters":
                ax.set_yscale('log')
            ax.legend(title=power_col)
    
        # Hide empty axes if n_metrics not filling grid
        for j in range(len(metrics), nrows * ncols):
            axes.flat[j].set_visible(False)
    
        plt.tight_layout()
    
        # Save or show
        if output_folder:
            os.makedirs(output_folder, exist_ok=True)
            filepath = os.path.join(output_folder, filename)
            plt.savefig(filepath, dpi=300, bbox_inches="tight")
            plt.close(fig)
            print(f"Plot saved to {filepath}")
            plt.show()
        else:
            plt.show()
            
    ## End helper functions--------------------------------------------------------------------

    ## Loop for clustering experiment
    results = pd.DataFrame(columns=["pcc_cutoff", "clustering_power", "mean_gene_margin", 
                                    "fraction_margin_positive", "modularity", "n_clusters", "n_genes_included"])

    if permuted_pcc_999 != None:
        pcc_cutoffs.append(permuted_pcc_999)

    # Run network construction and clustering for each hyperparam combination
    print("Running network clustering with parameters:")
    for i, pcc_cutoff in enumerate(pcc_cutoffs):
        print(f"PCC hard cutoff: {pcc_cutoff}, clustering soft power: ", end="")
    
        for j, clustering_power in enumerate(clustering_powers):
            print(f"{clustering_power}, ", end="")
    
            edge_list, node_label_df = make_spatial_network(pearsonR_mat,
                                                    gene_names=gene_names,
                                                    pcc_cutoff=pcc_cutoff,
                                                    clustering_power=clustering_power,
                                                    gene_labels_list=None,
                                                    gene_labels_names=None,
                                                    trials=infomap_clustering_trials,
                                                    output_folder=None)
            
            # Filter node_label_df and pearsonR_mat for non-singleton genes
            node_label_modules_df = node_label_df[
                node_label_df.groupby('module_label')['module_label'].transform('count') >= 2
            ]
            
            # Filter pearsonR mat and reorder according to node_label_modules_df
            gene_to_index = {g: i for i, g in enumerate(gene_names)}
            order_idx = [gene_to_index[g] for g in node_label_modules_df['name'] if g in gene_to_index]
            pearsonR_mat_filtered = pearsonR_mat[np.ix_(order_idx, order_idx)]

            # Calculate clustering metrics
            margin_scores, positive_margin, frac_positive, modularity, n_nodes, n_clusters = \
            get_clustering_stats(pearsonR_mat_filtered, node_label_modules_df, pcc_cutoff, clustering_power, min_genes_for_module)
    
            # Update results data frame
            results.loc[len(results)] = [pcc_cutoff, clustering_power, np.mean(margin_scores), 
                                         frac_positive, modularity, n_clusters, n_nodes]
        print("")

    ## Generate plot and save results
    if output_folder:
        plot_cutoff_grid(results, output_folder=output_folder, filename="clustering_plots.png")
        results.to_csv(f"{output_folder}/clustering_expt.csv", index=False)
    else:
        plot_cutoff_grid(results, output_folder=None, filename="clustering_plots.png")



def select_smoothing_degree(adata, micron_to_unit_conversion_factor, output_folder=None):
    """
    Evaluates a range of Gaussian smoothing standard deviations by measuring the spatial 
    frequency content of smoothed gene expression images, helping identify a suitable 
    gaussian_sd for run_parallelized_smoothing.

    Parameters:
    -----------
    adata : AnnData
        The input unsmoothed, library-size-normalized AnnData object with spatial coordinates 
        stored in adata.obsm['spatial'].

    micron_to_unit_conversion_factor : float
        Conversion factor from microns to the coordinate units used in adata.obsm['spatial'].
        This is platform-specific (e.g., 1.0 for Xenium, 2.0 for Stereo-seq).
        See docs/guides/micron_to_unit_conversion_table.md for platform-specific values.

    output_folder : str, optional
        Directory where output diagnostic plots will be saved.
        If None, plots are displayed but not saved. Default is None.

    Returns:
    --------
    None
        Displays diagnostic plots of spectral smoothness metrics (spectral_centroid, 
        spectral_slope, frac_high_freq) across Gaussian std. dev. values in microns 
        (5, 10, 20, 30, 40, 50, 60, 80, 100) plus an unsmoothed baseline (0).
        A good gaussian_sd is one where spectral metrics plateau, indicating diminishing 
        returns from further smoothing.
    """

    ## Helper functions--------------------------------------------------------------------

    def subsample_adata(adata, stride, target_n_output_spots=100000, n_genes=1000, min_UMIs_gene_filter=100):
        '''
        Params:
        adata: unsmoothed, normalized adata
        stride: stride for smoothing, in ST platform units
        target_n_output_spots: number of target spots in post-smoothing adata
        n_genes: number of genes to keep for smoothing (random sampling)
        '''

        def densest_square_region(adata, region_w):
            x = adata.obsm['spatial'][:, 0]
            y = adata.obsm['spatial'][:, 1]
        
            x_min, x_max = x.min(), x.max()
            y_min, y_max = y.min(), y.max()
        
            # clip region_w to the smaller dimension to ensure square
            region_w = min(region_w, x_max - x_min, y_max - y_min)
        
            n_x_bins = int(np.floor((x_max - x_min) / region_w))
            n_y_bins = int(np.floor((y_max - y_min) / region_w))
        
            # ensure we take the same number of bins along both axes
            n_bins = min(n_x_bins, n_y_bins)
        
            x_bins = x_min + np.arange(n_bins + 1) * region_w
            y_bins = y_min + np.arange(n_bins + 1) * region_w
        
            if len(x_bins) < 2 or len(y_bins) < 2:
                raise ValueError("region_w too large to fit inside spatial range.")
        
            hist, x_edges, y_edges = np.histogram2d(x, y, bins=[x_bins, y_bins])
        
            i_max, j_max = np.unravel_index(np.argmax(hist), hist.shape)
        
            x_center = (x_edges[i_max] + x_edges[i_max + 1]) / 2
            y_center = (y_edges[j_max] + y_edges[j_max + 1]) / 2
        
            spot_mask = (
                (x >= x_edges[i_max]) & (x < x_edges[i_max + 1]) &
                (y >= y_edges[j_max]) & (y < y_edges[j_max + 1])
            )
        
            return spot_mask, x_center, y_center
            

        # select a dense square region of tissue
        region_w = (target_n_output_spots ** 0.5) * stride
        spot_mask, x_center, y_center = densest_square_region(adata, region_w)
        adata._inplace_subset_obs(spot_mask)

        # update gene cols
        col_sums = np.sum(adata.X, axis=0)
        adata._inplace_subset_var(col_sums > min_UMIs_gene_filter)

        # subsample to n_genes many genes
        if adata.n_vars < n_genes:
            return adata  
        mask = np.zeros(adata.n_vars, dtype=bool)
        mask[np.random.choice(adata.n_vars, 1000, replace=False)] = True
        adata._inplace_subset_var(mask)
        
        return adata

    def square_bin(adata, stride):
        # Extract coordinates
        coords = adata.obsm['spatial']
        x, y = coords[:, 0], coords[:, 1]
        
        # Compute bin indices for each spot
        x_bins = np.floor(x / stride).astype(int)
        y_bins = np.floor(y / stride).astype(int)
        
        # Make a bin ID string (or tuple)
        bin_ids = [f"{i}_{j}" for i, j in zip(x_bins, y_bins)]
        
        # Create a DataFrame for grouping
        df = pd.DataFrame(adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X)
        df['bin_id'] = bin_ids
        
        # Group by bin and sum counts
        binned_expr = df.groupby('bin_id').sum(numeric_only=True)
        
        # Extract new spatial coordinates (bin centers)
        unique_bins = np.array([list(map(int, b.split('_'))) for b in binned_expr.index])
        bin_centers = unique_bins * stride + stride / 2.0  # center of the bin
    
        # Construct new AnnData
        binned_adata = ad.AnnData(
            X=binned_expr.values,
            var=adata.var.copy(),
            obs=pd.DataFrame(index=binned_expr.index)
        )
        binned_adata.obsm['spatial'] = bin_centers
    
        return binned_adata
        
    def add_pixel_coords(adata):
        
        # get sorted unique values for x and y
        x_unique = np.unique(adata.obsm['spatial'][:, 0])
        y_unique = np.unique(adata.obsm['spatial'][:, 1])
        
        # build mapping dicts
        x_coord_dict = {val: i for i, val in enumerate(x_unique)}
        y_coord_dict = {val: i for i, val in enumerate(y_unique)}
        
        # map spatial coords to integer pixel indices
        pixel_indices = np.column_stack([
            [x_coord_dict[x] for x in adata.obsm['spatial'][:, 0]],
            [y_coord_dict[y] for y in adata.obsm['spatial'][:, 1]]
        ])
        
        # store in obsm
        adata.obsm['pixel_indices'] = pixel_indices
    
        return
    
    
    def generate_img_stack(adata):
    
        x = adata.obsm['pixel_indices'][:,0]
        y = adata.obsm['pixel_indices'][:,1]
        
        width = len(np.unique(x))
        height = len(np.unique(y))
        
        gene_img_stack = np.zeros((len(adata.var_names), height, width), dtype=np.float32)
        
        for i, gene in enumerate(adata.var_names):
            
            col_idx = adata.var_names.get_loc(gene)
            gene_vals = adata.X[:,col_idx]
        
            gene_img_stack[i, y, x] = gene_vals
    
        return gene_img_stack
    
    
    def generate_tissue_region_mask(adata):
    
        total_gene_vals = np.sum(adata.X, axis=1)
    
        x = adata.obsm['pixel_indices'][:,0]
        y = adata.obsm['pixel_indices'][:,1]
        
        width = len(np.unique(x))
        height = len(np.unique(y))
        
        total_UMI_img = np.zeros((height, width), dtype=total_gene_vals.dtype)
        total_UMI_img[y, x] = total_gene_vals
        
        tissue_region_mask = total_UMI_img > 0
    
        return tissue_region_mask

        
    def compute_channel_smoothness_metrics(img_stack,
                                           tissue_region_mask,
                                           microns_gaussian_sd,
                                           microns_per_pixel,
                                           window=True,
                                           normalize=True,
                                           radial_max_frac=0.5):
        """
        img_stack: array of shape (C, H, W) or (H, W) for single channel.
        Returns: pandas.DataFrame with metrics per channel.
        """
        # ensure shape (C, H, W)
        single = False
        if img_stack.ndim == 2:
            img_stack = img_stack[None, ...]
            single = True
        C, H, W = img_stack.shape
    
        # optional windowing (2D Hann)
        if window:
            wx = np.hanning(W)
            wy = np.hanning(H)
            win2d = np.outer(wy, wx)
        else:
            win2d = np.ones((H, W), dtype=float)
    
        # precompute radius bins for radial profile (centered)
        cy, cx = H // 2, W // 2
        y, x = np.indices((H, W))
        r = np.hypot(x - cx, y - cy)
        r_int = r.astype(int)
        max_r = int(np.floor(min(cx, cy)))  # meaningful radial range
        bins = np.arange(0, max_r + 1)
    
        # precompute masks for spectral slope fitting (freqs to consider)
        # avoid DC (bin 0) and very high bins near Nyquist
        radial_idxs_for_slope = np.arange(1, int(radial_max_frac * len(bins)) + 1)
    
        rows = []
        eps = 1e-12
    
        for ch in range(C):
            img = img_stack[ch].astype(float)
    
            mean_nonzero = img[tissue_region_mask].mean()
            img[tissue_region_mask] -= mean_nonzero
    
            # optional normalize by RMS or STD to compare across channels
            if normalize:
                denom = np.sqrt(np.nanmean(img**2)) + eps  # RMS
                img = img / denom
    
            # apply window
            img_win = img * win2d
    
            # spectral-domain metrics
            F = np.fft.fft2(img_win)
            F = np.fft.fftshift(F)
            P = np.abs(F)**2
    
            # remove DC (center) explicitly
            P[cy, cx] = 0.0
    
            # radial profile using bincount
            tbin = np.bincount(r_int.ravel(), weights=P.ravel(), minlength=len(bins))
            nr = np.bincount(r_int.ravel(), minlength=len(bins))
            with np.errstate(divide='ignore', invalid='ignore'):
                radial_mean = tbin / np.maximum(nr, 1)  # shape (len(bins),)
            freqs = np.arange(len(radial_mean)) / (W * microns_per_pixel)
            
            # fraction of power spectrum with frequency above 0.025 cycles/micron (noise frequency levels)
            numer_high_frac = radial_mean[freqs > 0.025].sum()
            denom_high_frac = radial_mean.sum() + eps
            high_frac = numer_high_frac / denom_high_frac
    
            # spectral centroid (weighted by radial_mean)
            denom_centroid = radial_mean.sum() + eps
            spectral_centroid = (freqs * radial_mean).sum() / denom_centroid
    
            # spectral slope: fit log-log to radial_mean over selected bins (avoid zeros)
            fit_idx = radial_idxs_for_slope
            vals = radial_mean[fit_idx]
            xfit = freqs[fit_idx]
            mask = (vals > 0) & (xfit > 0)
            if mask.sum() >= 2:
                slope, intercept = np.polyfit(np.log(xfit[mask]), np.log(vals[mask]), 1)
            else:
                slope = np.nan
    
            rows.append({
                "gaussian_sd": microns_gaussian_sd,
                "frac_high_freq": high_frac,
                "spectral_centroid": spectral_centroid,
                "spectral_slope": slope
            })
    
        df = pd.DataFrame(rows)
        if single:
            return df.iloc[0].to_dict()
        return df
        

    def plot_results(df, output_folder=None, filename="smoothness_expt_plots_split.png"):
        metrics = ["spectral_centroid", "spectral_slope", "frac_high_freq"]
        n_metrics = len(metrics)
    
        fig, axes = plt.subplots(1, n_metrics, figsize=(6 * n_metrics, 4))
    
        for ax, metric in zip(axes, metrics):

            if metric != "spectral_slope":
                # Split data: gaussian_sd==0 vs gaussian_sd!=0
                df_zero = df[df['gaussian_sd'] == 0]
                df_nonzero = df[df['gaussian_sd'] != 0]
        
                sns.violinplot(data=df_nonzero, x="gaussian_sd", y=metric, inner="box", cut=0, ax=ax, color="skyblue", edgecolor="black")
        
                if not df_zero.empty:
                    max_nonzero = df_nonzero[metric].max()
                    max_zero = df_zero[metric].max()
                    # create secondary y-axis for the Gaussian_SD = 0 rows
                    ax2 = ax.twinx()
                    sns.violinplot(data=df_zero, x="gaussian_sd", y=metric, inner="box", cut=0, ax=ax2, color="salmon", edgecolor="black")
                    # scale the secondary axis to show only high values
                    ax2.set_ylim(max_nonzero * 1.5, max_zero * 1.1)
                    ax2.set_ylabel(f"{metric} (Gaussian std. dev. = 0)")
                    ax2.grid(False)
                
                ax.set_title(metric)
                ax.set_xlabel("Gaussian std. dev. (microns)")
                ax.set_ylabel(metric)
            
            else: # metric == 'spectral_slope'
                sns.violinplot(data=df, x="gaussian_sd", y=metric, inner="box", cut=0, ax=ax, color="skyblue", edgecolor="black") 
                ax.set_title(metric)
                ax.set_xlabel("Gaussian std. dev. (microns)")
                ax.set_ylabel(metric)
            
        plt.tight_layout()
    
        if output_folder:
            import os
            os.makedirs(output_folder, exist_ok=True)
            save_path = os.path.join(output_folder, filename)
            plt.savefig(save_path, dpi=150)
            print(f"Saved to: {save_path}")
    
        plt.show()
        plt.close()



    ## End helper functions----------------------------------------------------------------

    
    gaussian_SDs_microns = [5,10,20,30,40,50,60,80,100, 0] # testing values in microns
    gaussian_SDs = list(np.array(gaussian_SDs_microns) * micron_to_unit_conversion_factor)

    ## Subset adata for efficiency
    s_adata = subsample_adata(adata.copy(), stride=gaussian_SDs[0])
    
    fft_df_list = []
    
    ## Run smoothing expt and collect smoothness measure across gaussian_sd inputs
    for i, gaussian_sd in enumerate(gaussian_SDs):

        if gaussian_sd > 0:
            # Run Smoothing
            print(f"Smoothing with Gaussian std. dev. = {gaussian_SDs_microns[i]} microns...")
            sm_s_adata = run_gridbased_smoothing(
                ADATA=s_adata,
                STRIDE=gaussian_SDs[0],
                GRID_FITTING_DIST=-1, # not used here
                MIN_BEADS_UNDER_KERNEL=1, # FIX THIS
                GAUSSIAN_SD=gaussian_sd, 
                KERNEL_RADIUS=3*gaussian_sd,
                NUM_PROCESSES=10,
                NUM_DATA_SPLITS=2,
                IS_FFT_EXPT=True
            )
        else: # gaussian_sd == 0
            sm_s_adata = square_bin(adata=s_adata, stride=gaussian_SDs[0])
        
        # Processing prior to FFT
        add_pixel_coords(sm_s_adata)
        gene_img_stack = generate_img_stack(sm_s_adata)
        tissue_region_mask = generate_tissue_region_mask(sm_s_adata)
    
        # Run FFT to get smoothness measurements
        df = compute_channel_smoothness_metrics(
            gene_img_stack,
            tissue_region_mask,
            microns_per_pixel=gaussian_SDs_microns[0],
            microns_gaussian_sd=gaussian_SDs_microns[i],
            window=True,
            normalize=True,
            radial_max_frac=0.5
        )

        fft_df_list.append(df)

    ## Output experiment results plots / table
    final_df = pd.concat(fft_df_list, axis=0, ignore_index=True)
    
    plot_results(final_df, output_folder=output_folder)
    