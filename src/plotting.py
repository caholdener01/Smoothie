# Copyright (c) 2025 Chase Holdener
# Licensed under the MIT License. See LICENSE file for details.

import os
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

### Plotting functions for single dataset pipeline

def plot_gene(sm_adata, gene_name, output_folder, figsize=(1,1), 
              spot_size=25, fontsize=6, fontfamily='sans-serif', cmap=None,
              dpi=300, save_plot=True):
    """
    Plots the spatial expression of a specified gene in a single spatial transcriptomic dataset.

    Parameters
    ----------
    sm_adata : AnnData 
        The AnnData object containing spatial transcriptomic data.
    
    gene_name : str
        Name of the gene to visualize.
    
    output_folder : str
        Directory path where the output plot will be saved.
    
    figsize : tuple, optional
        Figure size in inches (width, height). Default is (1,1).
    
    spot_size : int, optional
        Size of spatial transcriptomic spots in the plot. Default is 25.
    
    fontsize : int, optional
        Font size for the plot title. Default is 6.
    
    fontfamily : str, optional
        Font family used for the plot title. Default is 'sans-serif'.
    
    cmap : matplotlib.colors.Colormap, optional
        Colormap to use for visualizing spatial expression patterns.
        If None, defaults to a custom gray-red-black colormap.
    
    dpi : int, optional
        Resolution of the saved plot in dots per inch (DPI). Default is 1200.
    
    save_plot : bool, optional
        If True, saves the generated plot as a pdf file in the specified output folder. Default is True.
    
    Returns
    -------
    None
        Displays the plot and optionally saves it to the output folder.
    """
    sc.settings.set_figure_params(dpi=150, facecolor="white")
    
    if cmap is None:
        colors = ['#eaeaea', '#bb000b', 'black']
        cmap = LinearSegmentedColormap.from_list("custom_gray_red_black", colors)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    sc.pl.spatial(
        sm_adata,
        color=[gene_name],
        cmap=cmap,
        vmax=None,
        colorbar_loc=None,
        title=gene_name,
        spot_size=spot_size,
        show=False,
        ax=ax
    )
    ax.set_axis_off()
    ax.set_title(f'{gene_name}', fontsize=fontsize, fontfamily=fontfamily, y=0.95)
    plt.tight_layout(pad=0) 

    if save_plot:
        pdf_filename = f"{output_folder}/{gene_name}.pdf"
        plt.savefig(pdf_filename, format='pdf', bbox_inches='tight', dpi=dpi, pad_inches=0)
    plt.show()


def plot_modules(sm_adata, node_label_df, output_folder, plots_per_row, min_genes=3, 
                 figsize=None, spot_size=25, fontsize=6, fontfamily='sans-serif', 
                 cmap=None, dpi=300, save_plot=True):
    """
    Plots module scores using sc.pl.spatial, arranging plots in rows with a fixed number of columns.

    Parameters
    ----------
    sm_adata : AnnData
        The AnnData object containing spatial transcriptomic data.

    node_label_df : pd.DataFrame
        DataFrame containing gene module assignments with at least two columns:
        'community_label' (gene module ID) and 'name' (gene name).
    
    output_folder : str
        Directory path where the output plots will be saved.

    plots_per_row : int
        Number of plots per row in each output pdf.

    min_genes : int, optional
        Minimum number of genes required for a module to be plotted. Default is 3.

    figsize : tuple, optional
        Tuple specifying the figure size in inches (width, height). Default is (plots_per_row, 1).

    spot_size : int, optional
        Size of spatial transcriptomic spots in the plots. Default is 25.

    fontsize : int, optional
        Font size for subplot titles. Default is 6.

    fontfamily : str, optional
        Font family used for the plot titles. Default is 'sans-serif'.

    cmap : matplotlib.colors.Colormap, optional
        Colormap to use for visualizing spatial expression patterns.
        If None, defaults to a custom gray-red-black colormap.

    dpi : int, optional
        Resolution of saved plot images in dots per inch (DPI). Default is 1200.

    save_plot : bool, optional
        If True, saves the generated plots as pdf files in the specified output folder. Default is True.

    Returns
    -------
    None
        Displays the plots and optionally saves them to the output folder.
    """
    sc.settings.set_figure_params(dpi=150, facecolor="white")
    
    if figsize is None:
        figsize = (plots_per_row, 1)
        
    if cmap is None:
        colors = ['#eaeaea', '#bb000b', 'black']
        cmap = LinearSegmentedColormap.from_list("custom_gray_red_black", colors)
    
    # Extract and filter gene modules
    gene_modules = node_label_df.groupby('community_label')['name'].apply(list).tolist()
    gene_modules = [gene_list for gene_list in gene_modules if len(gene_list) >= min_genes]
    
    # Determine number of rows needed
    total_modules = len(gene_modules)
    num_rows = int(np.ceil(total_modules / plots_per_row))

    for row in range(num_rows):
        start_idx = row * plots_per_row
        end_idx = min(start_idx + plots_per_row, total_modules)
        modules_to_plot = gene_modules[start_idx:end_idx]
        
        fig, axes = plt.subplots(nrows=1, ncols=len(modules_to_plot), figsize=figsize)
        
        if len(modules_to_plot) == 1:
            axes = [axes]

        counter=1
        for ax, module_genes in zip(axes, modules_to_plot):

            # Compute module score (sum of 0-to-1 rescaled gene expression from the module)
            temp_gene_idx = [gene in module_genes for gene in sm_adata.var_names]
            temp_sm_adata = sm_adata.X[:, temp_gene_idx]
            col_maxs = np.max(temp_sm_adata, axis=0)
            temp_norm_sm_adata = temp_sm_adata / col_maxs
            combined_cluster_surface = np.sum(temp_norm_sm_adata, axis=1)
            sm_adata.obs['module_score'] = combined_cluster_surface

            sc.pl.spatial(
                sm_adata,
                color=['module_score'],
                cmap=cmap,
                vmax=None,
                title='',
                colorbar_loc=None,
                spot_size=spot_size,
                show=False,
                ax=ax
            )
            ax.set_title(f"Module {start_idx + counter}", fontsize=fontsize, fontfamily=fontfamily, y=0.95)
            ax.set_axis_off()
            
            counter += 1
        
        plt.subplots_adjust(wspace=0, hspace=0)

        if save_plot:
            os.makedirs(output_folder, exist_ok=True)
            pdf_filename = f"{output_folder}/module_scores_row_{row+1}.pdf"
            plt.savefig(pdf_filename, format='pdf', bbox_inches='tight', dpi=dpi, pad_inches=0)
        
        plt.show()




### Plotting functions for multi dataset pipeline

def plot_gene_multisample(sm_adata_set, adata_set_names, gene_name, output_folder, shared_scaling=False, 
                          figsize=None, spot_size=25, fontsize=6, fontfamily='sans-serif', cmap=None,
                          dpi=300, save_plot=True):
    """
    Plots the spatial expression of a specified gene across multiple spatial transcriptomic datasets,
    with an option for shared scaling across datasets.

    Parameters
    ----------
    sm_adata_set : list of AnnData
        List of AnnData objects representing spatial transcriptomic datasets.
        
    adata_set_names : list of str
        Names of the datasets for labeling subplots.
    
    gene_name : str
        Name of the gene to be visualized.
    
    output_folder : str
        Directory where the output plots will be saved.

    shared_scaling : bool, optional
        If True, normalizes gene expression using the maximum value of the gene across all datasets.
        If False, normalizes independently within each dataset. Default is False.
        
    figsize : tuple, optional
        Tuple specifying the figure size in inches (width, height). If None, defaults to (NCOLS, 1).
    
    spot_size : int, optional
        Size of spatial transcriptomic spots in the plots. Default is 25.
    
    fontsize : int, optional
        Font size for subplot titles. Default is 6.
    
    fontfamily : str, optional
        Font family used for text annotations. Default is 'sans-serif'.
    
    cmap : matplotlib.colors.Colormap, optional
        Colormap to use for visualizing spatial expression patterns. 
        If None, defaults to a custom gray-red-black colormap.
    
    dpi : int, optional
        Resolution of saved plot images in dots per inch (DPI). Default is 1200.
    
    save_plot : bool, optional
        If True, saves the generated plot as a pdf file in the specified output folder. Default is True.
    
    Returns
    -------
    None
        Displays the plots and optionally saves them.
    """
    sc.settings.set_figure_params(dpi=150, facecolor="white")
    
    if figsize is None:
        figsize = (len(sm_adata_set), 1)

    if cmap is None:
        colors = ['#eaeaea', '#bb000b', 'black']
        cmap = LinearSegmentedColormap.from_list("custom_gray_red_black", colors)

    # Compute spatial scaling factors
    max_x_range, max_y_range = 0.0, 0.0
    center_points = []

    for adata in sm_adata_set:
        x_coords, y_coords = adata.obsm['spatial'][:, 0], adata.obsm['spatial'][:, 1]
        x_min, x_max = x_coords.min(), x_coords.max()
        y_min, y_max = y_coords.min(), y_coords.max()
        
        x_range, y_range = x_max - x_min, y_max - y_min
        max_x_range, max_y_range = max(max_x_range, x_range), max(max_y_range, y_range)

        center_x, center_y = x_min + 0.5 * x_range, y_min + 0.5 * y_range
        center_points.append([center_x, center_y])

    # Determine shared max expression value if using shared scaling
    global_max = None
    if shared_scaling:
        global_max = 0
        for adata in sm_adata_set:
            if gene_name in adata.var_names:
                gene_idx = adata.var_names.get_loc(gene_name)
                global_max = max(global_max, adata.X[:, gene_idx].max())

    # Create subplots
    fig, axes = plt.subplots(nrows=1, ncols=len(sm_adata_set), figsize=figsize)

    for i, (adata, ax) in enumerate(zip(sm_adata_set, axes)):
        if gene_name in adata.var_names:
            gene_idx = adata.var_names.get_loc(gene_name)
            vmax_value = global_max if shared_scaling else None

            sc.pl.spatial(adata, 
                          color=[gene_name], 
                          cmap=cmap, 
                          title='', 
                          colorbar_loc=None, 
                          spot_size=spot_size, 
                          show=False, 
                          ax=ax, 
                          vmax=vmax_value)
        else:
            # If gene is missing, plot an empty (constant color) image
            sc.pl.spatial(adata, 
                          color=[adata.var_names[0]], 
                          cmap=cmap, 
                          vmin=1e7, 
                          vmax=1e7+1, 
                          title='', 
                          colorbar_loc=None, 
                          spot_size=spot_size, 
                          show=False, 
                          ax=ax)

        ax.set_title(f'{adata_set_names[i]}', fontsize=fontsize, fontfamily=fontfamily, y=0.95)

        # Set consistent spatial scaling
        x_min, x_max = center_points[i][0] - 0.5 * max_x_range, center_points[i][0] + 0.5 * max_x_range
        y_min, y_max = center_points[i][1] - 0.5 * max_y_range, center_points[i][1] + 0.5 * max_y_range
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_max, y_min)
        ax.set_axis_off()

    plt.subplots_adjust(wspace=0, hspace=0)

    if save_plot:
        os.makedirs(output_folder, exist_ok=True)
        scaling_tag = "_shared_scale" if shared_scaling else "_independent_scale"
        pdf_filename = f'{output_folder}/{gene_name}{scaling_tag}.pdf'
        plt.savefig(pdf_filename, format='pdf', bbox_inches='tight', dpi=dpi, pad_inches=0)

    plt.show()



def plot_modules_multisample(sm_adata_set, adata_set_names, node_label_df, output_folder, 
                             shared_scaling=False, min_genes=3, figsize=None, spot_size=25,
                             fontsize=6, fontfamily='sans-serif', cmap=None, dpi=300, save_plot=True):
    """
    Plots multi-dataset gene modules with an option for shared scaling across datasets.
    
    Parameters:    
    -----------
    sm_adata_set (list of AnnData): 
        List of AnnData objects representing spatial transcriptomic datasets.
        
    adata_set_names (list of str): 
        List of dataset names for labeling subplots.
    
    node_label_df (pd.DataFrame): 
        DataFrame containing gene module assignments with at least two columns: 
        'community_label' (gene module ID) and 'name' (gene name).
    
    output_folder (str): 
        Directory path where output plots will be saved.
        
    shared_scaling (bool, optional): 
        If True, normalizes gene expression using the maximum value of each gene across all datasets.
        If False, normalizes each gene independently within each dataset. Default is False.
        
    min_genes (int, optional): 
        Minimum number of genes required for a module to be plotted. Default is 3.
        
    figsize (tuple, optional): 
        Tuple specifying the figure size in inches (width, height). If None, defaults to (NCOLS, 1). 
    
    spot_size (int, optional): 
        Size of spatial transcriptomic spots in the plots. Default is 25.
    
    fontsize (int, optional): 
        Font size for subplot titles. Default is 6.
    
    fontfamily (str, optional): 
        Font family used for text annotations. Default is 'sans-serif'.
    
    cmap (matplotlib.colors.Colormap, optional): 
        Colormap to use for visualizing spatial expression patterns. 
        If None, defaults to a custom gray-red-black colormap.
    
    dpi (int, optional): 
        Resolution of saved plot images in dots per inch (DPI). Default is 1200.
    
    save_plot (bool, optional): 
        If True, saves the generated plots as pdf files in the specified output folder. Default is True.
    
    Returns:
    --------
    None: 
        Displays the plots and optionally saves them to the output folder.
    """
    sc.settings.set_figure_params(dpi=150, facecolor="white")
    
    if figsize == None:
        figsize = (len(sm_adata_set), 1)
    if cmap == None:
        colors = ['#eaeaea' , '#bb000b', 'black']
        custom_cmap = LinearSegmentedColormap.from_list("custom_gray_red_black", colors)
        cmap = custom_cmap
        
    # Extract gene modules
    gene_modules = node_label_df.groupby('community_label')['name'].apply(list).tolist()
    gene_modules = [gene_list for gene_list in gene_modules if len(gene_list) >= min_genes]
    
    # Determine spatial ranges for consistent tissue scaling
    max_x_range, max_y_range = 0.0, 0.0
    center_points = []
    for i in range(len(sm_adata_set)):
        x_max, x_min = sm_adata_set[i].obsm['spatial'][:, 0].max(), sm_adata_set[i].obsm['spatial'][:, 0].min()
        y_max, y_min = sm_adata_set[i].obsm['spatial'][:, 1].max(), sm_adata_set[i].obsm['spatial'][:, 1].min()
        x_range, y_range = x_max - x_min, y_max - y_min
        max_x_range, max_y_range = max(max_x_range, x_range), max(max_y_range, y_range)
        center_points.append([x_min + (0.5 * x_range), y_min + (0.5 * y_range)])

    for i, temp_geneset in enumerate(gene_modules):

        if shared_scaling:
            # Shared scaling: Find max expression values across all datasets
            all_adata_col_maxs = np.zeros(len(temp_geneset))
            for d, temp_adata in enumerate(sm_adata_set):
                temp_geneset_idx = [temp_adata.var_names.get_loc(gene) if gene in temp_adata.var_names else None for gene in temp_geneset]
                temp_sm_adata = temp_adata.X  # Get expression matrix
                for j, idx in enumerate(temp_geneset_idx):
                    if idx is not None:
                        all_adata_col_maxs[j] = max(all_adata_col_maxs[j], np.max(temp_sm_adata[:, idx]))

            # Normalize each dataset using shared max values
            for d, temp_adata in enumerate(sm_adata_set):
                temp_geneset_idx = [temp_adata.var_names.get_loc(gene) if gene in temp_adata.var_names else None for gene in temp_geneset]
                temp_sm_adata = temp_adata.X
                normalized_sums = np.zeros(temp_adata.n_obs)
                for j, idx in enumerate(temp_geneset_idx):
                    if idx is not None:
                        normalized_col = temp_sm_adata[:, idx] / all_adata_col_maxs[j]
                        normalized_sums += np.nan_to_num(normalized_col)
                temp_adata.obs['normalized_sum_surface'] = normalized_sums

            # Determine vmax across all datasets
            vmax_val = max(np.max(temp_adata.obs['normalized_sum_surface']) for temp_adata in sm_adata_set)
        
        # Create figure for module
        fig, axes = plt.subplots(nrows=1, ncols=len(sm_adata_set), figsize=figsize)
        for d, temp_adata in enumerate(sm_adata_set):
            ax = axes[d]

            if shared_scaling:
                color_key = 'normalized_sum_surface'
            else:
                # Independent scaling: Normalize each dataset separately
                temp_gene_idx = [gene in temp_geneset for gene in temp_adata.var_names]
                temp_sm_adata = temp_adata.X[:, temp_gene_idx]
                col_maxs = np.max(temp_sm_adata, axis=0)
                temp_norm_sm_adata = temp_sm_adata / col_maxs
                combined_cluster_surface = np.sum(temp_norm_sm_adata, axis=1)
                temp_adata.obs['temp_combined_cluster_surface'] = combined_cluster_surface
                color_key = 'temp_combined_cluster_surface'
                vmax_val = None

            sc.pl.spatial(temp_adata, 
                          color=[color_key],
                          cmap=custom_cmap, 
                          vmax=vmax_val, 
                          title='',
                          colorbar_loc=None,
                          spot_size=spot_size,
                          show=False,
                          ax=ax)
            
            ax.set_axis_off()
            ax.set_title(f'{adata_set_names[d]}', fontsize=fontsize, fontfamily=fontfamily, y=0.95)
            x_min, x_max = center_points[d][0] - (0.5 * max_x_range), center_points[d][0] + (0.5 * max_x_range)
            y_min, y_max = center_points[d][1] - (0.5 * max_y_range), center_points[d][1] + (0.5 * max_y_range)
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_max, y_min)

        plt.subplots_adjust(wspace=0, hspace=0)
        if save_plot:
            os.makedirs(output_folder, exist_ok=True)
            scaling_tag = "_shared_scale" if shared_scaling else "_independent_scale"
            pdf_filename = f'{output_folder}/Module_{i+1}{scaling_tag}.pdf'
            plt.savefig(pdf_filename, format='pdf', bbox_inches='tight', dpi=dpi, pad_inches=0)
        plt.show()
