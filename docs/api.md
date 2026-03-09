# Smoothie API Reference

---

## gaussian_smoothing.py

### `run_parallelized_smoothing`

```python
run_parallelized_smoothing(adata, grid_based_or_not, gaussian_sd, min_spots_under_gaussian,
                           stride=None, grid_fitting_dist=None, num_processes=10, num_data_splits=None)
```

Performs parallelized Gaussian smoothing on spatial transcriptomics data.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `adata` | AnnData | The input dataset. |
| `grid_based_or_not` | bool | `True` = grid-based smoothing (smooth only at hexagonal grid points); `False` = in-place smoothing (smooth at every spatial location). In-place is effective for cell-sized or larger resolution (10–50 µm). Grid-based is effective for subcellular to cell-sized resolution (0.5–10 µm). |
| `gaussian_sd` | float | Standard deviation for the Gaussian kernel in coordinate units. Use `target_microns × micron_to_unit_conversion`. Target of 20–40 µm is generally ideal. See `docs/guides/micron_to_unit_conversion_table.md`. |
| `min_spots_under_gaussian` | int | Minimum number of data points within radius `3 × gaussian_sd` required for smoothing to occur at a location (typical range: 25–100). |
| `stride` | float, optional | Stride value for grid-based smoothing. Default is `1 × gaussian_sd`; `0.5 × gaussian_sd` allows a denser grid. |
| `grid_fitting_dist` | float, optional | Minimum distance a grid point must be from a spatial coordinate to be retained. Default is `0.25 × gaussian_sd`. |
| `num_processes` | int, optional | Number of parallel processes to use. Default is 10. |
| `num_data_splits` | int, optional | Number of data chunks for memory-efficient processing. If unspecified, automatically selected. |

**Returns:** `sm_adata` (AnnData) — The smoothed AnnData object.

---

## spatial_correlation.py

### `compute_correlation_matrix`

```python
compute_correlation_matrix(X)
```

Computes pairwise Pearson correlation coefficients between all genes.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `X` | np.ndarray | An `(N × G)` matrix where rows are spatial points and columns are genes. Should be smoothed gene expression data. |

**Returns:** `(pearsonR_mat, p_val_mat)` — Tuple of two `(G × G)` matrices: Pearson correlation coefficients and corresponding p-values. 

---

### `get_correlations_to_GOI`

```python
get_correlations_to_GOI(pearsonR_mat, gene_names, GOI, reverse_order=False, plot_histogram=True)
```

Retrieves and ranks the correlation of all genes with a specified gene of interest (GOI).

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `pearsonR_mat` | np.ndarray | A `(G × G)` Pearson correlation coefficient matrix. |
| `gene_names` | list, np.ndarray, or pd.Index | Gene names corresponding to matrix indices. |
| `GOI` | str | The gene of interest for which correlations are ranked. |
| `reverse_order` | bool, optional | If `True`, sorts correlations in ascending order. Default is `False` (descending). |
| `plot_histogram` | bool, optional | If `True`, plots a histogram of correlation values to GOI. Default is `True`. |

**Returns:** `(G × 2)` np.ndarray — First column contains gene names, second column contains correlation values with the GOI, sorted by correlation strength.

---

## network_analysis.py

### `make_spatial_network`

```python
make_spatial_network(pearsonR_mat, gene_names, pcc_cutoff, clustering_power,
                     output_folder=None, save_file_prefix="", gene_labels_list=None,
                     gene_labels_names=None, trials=20, random_seed=0)
```

Generates a spatial co-expression network with hard thresholding (`pcc_cutoff`) and soft power transformation (`clustering_power`), then performs Infomap clustering and calculates network metrics.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `pearsonR_mat` | np.ndarray | Square matrix of Pearson correlation coefficients between genes. |
| `gene_names` | list or pd.Index | Gene names corresponding to rows/columns of the correlation matrix. |
| `pcc_cutoff` | float | Threshold for the Pearson correlation coefficient; only correlations above this value are retained. |
| `clustering_power` | float | Soft power that controls rescaling of PCC values. Higher values favor more gene modules. |
| `output_folder` | str, optional | Directory where output files are saved. If `None`, no files are saved. Default is `None`. |
| `save_file_prefix` | str, optional | String prefix added to output filenames. |
| `gene_labels_list` | list of list/tuple/np.ndarray, optional | Gene set labels for visualization. Each item must have the same length as the number of genes. |
| `gene_labels_names` | list, optional | Names for the gene sets in `gene_labels_list`. |
| `trials` | int, optional | Number of Infomap clustering trials. Default is 20. |
| `random_seed` | int, optional | Random seed for the Infomap clustering algorithm. Default is 0. |

**Returns:** `(edge_list, node_label_df)` — Edge list in the format `[gene1, gene2, PCC, Rescaled_PCC]`, and a DataFrame with gene names, module labels, and network metrics.

---

### `make_geneset_spatial_network`

```python
make_geneset_spatial_network(pearsonR_mat, gene_names, node_label_df, gene_list,
                              low_pcc_cutoff, output_folder=None, intra_geneset_edges_only=False)
```

Given a clustered network, constructs a new network with a lower PCC cutoff applied to a provided gene set.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `pearsonR_mat` | np.ndarray | Pearson correlation coefficient matrix. |
| `gene_names` | list | Gene names corresponding to rows/columns of `pearsonR_mat`. |
| `node_label_df` | pd.DataFrame | DataFrame with gene names and module labels from `make_spatial_network`. |
| `gene_list` | list | Genes to retain in the subset network. |
| `low_pcc_cutoff` | float | Minimum PCC required to include an edge. Must be ≤ the `pcc_cutoff` used to generate `node_label_df`. |
| `output_folder` | str, optional | Directory where output files are saved. If `None`, no files are saved. Default is `None`. |
| `intra_geneset_edges_only` | bool, optional | If `True`, only includes edges where both nodes are in `gene_list`. Default is `False`. |

**Returns:** `(geneset_edge_list, geneset_node_label_df)` — Edge list in format `[gene1, gene2, PCC]`, and a DataFrame with updated module labels including column `weak_module_label` for gene_list genes that have a correlation above low_pcc_cutoff with another gene in the previously clustered network (node_label_df modules).

---

## multi_sample_integration.py

### `concatenate_smoothed_matrices`

```python
concatenate_smoothed_matrices(sm_adata_list)
```

Concatenates smoothed count matrices from multiple datasets, aligning gene columns across datasets.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `sm_adata_list` | list of AnnData | AnnData objects with smoothed count matrices. |

**Returns:** `(X_concat, gene_names)` — Concatenated `(N × G)` matrix with aligned genes, and the corresponding list of gene names.

---

### `run_second_order_correlation_analysis`

```python
run_second_order_correlation_analysis(sm_adata_list, pcc_cutoff, node_label_df=None,
                                      E_max=25, output_folder=None, seed=0)
```

Comprehensive function for second-order correlation gene embedding analysis across multiple samples. Workflow: (1) aligns correlation matrices, (2) creates gene embeddings, (3) computes stability metrics, and (4) adds stability values to `node_label_df`. Check this function's Returns section for information on output.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `sm_adata_list` | list of AnnData | Smoothed AnnData objects from different spatial transcriptomics samples. Each must have gene expression data in `.X` and gene names in `.var_names`. |
| `pcc_cutoff` | float | Minimum PCC threshold. Genes must exceed this value in at least one dataset to be included in embedding analysis. Typical values: 0.4–0.5. |
| `node_label_df` | pd.DataFrame, optional | Gene module assignments with columns `'name'` and `'module_label'`. If provided, embedding features are balanced across modules and stability metrics are added as a `'gene_stabilities'` column. Default is `None`. |
| `E_max` | int, optional | Maximum number of genes per module used for embedding features. Controls embedding dimensionality and balances representation across modules. Typical values: 15–50. Default is 25. |
| `output_folder` | str, optional | Path to folder where results are saved. If `None`, files are not saved. Default is `None`. |
| `seed` | int, optional | Random seed for reproducibility when downsampling embedding features. Default is 0. |

**Returns:** `dict` with keys:

- `'gene_embeddings_tensor'` (np.ndarray): Shape `(n_genes, n_features, n_datasets)`. Each gene's embedding is a vector of Fisher Z-transformed correlation values with other genes (features).
- `'robust_gene_names'` (np.ndarray): Shape `(n_genes)`. These are the genes that passed the `pcc_cutoff` correlation threshold in at least one dataset in sm_adata_list.
- `'gene_stabilities'` (np.ndarray): Shape `(n_datasets, n_datasets, n_genes)` Tensor of gene stabilities across dataset pairs. Values range from -1 to 1, where higher values indicate more stable/conserved co-expression patterns. A stability of -1 designates that a gene was missing in one of the datasets.
- `'aligned_pearsonR_tensor'` (np.ndarray): Shape `(n_all_genes, n_all_genes, n_datasets)`. Contains aligned Pearson correlation matrices across all datasets.
- `'node_label_df'` (pd.DataFrame or None): Input `node_label_df` with `'gene_stabilities'` column added, or `None`.

---

### `plot_dataset_pair_stabilities`

```python
plot_dataset_pair_stabilities(gene_stabilities, dataset_names, output_folder=None,
                               figsize=None, fontsize=5, dpi=300, file_format='png')
```

Visualizes the distribution of gene stability between all dataset pairs using violin plots. See `run_second_order_correlation_analysis` for more info.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `gene_stabilities` | np.ndarray | Stabilities tensor of shape `(n_datasets, n_datasets, n_genes)` from `run_second_order_correlation_analysis`. |
| `dataset_names` | list of str | Dataset names for axis labels. |
| `output_folder` | str, optional | Directory to save the plot. If `None`, does not save. |
| `figsize` | tuple, optional | Figure size in inches. If `None`, auto-calculated based on number of pairs. |
| `fontsize` | int, optional | Font size for labels and titles. Default is 5. |
| `dpi` | int, optional | Resolution for saved figure. Default is 300. |
| `file_format` | str, optional | `'png'` or `'pdf'`. Default is `'png'`. |

**Returns:** None.

---

### `plot_gene_stability_distribution`

```python
plot_gene_stability_distribution(gene_stabilities, robust_gene_names, node_label_df=None,
                                  output_folder=None, figsize=(3,2), fontsize=6,
                                  bins=50, dpi=300, file_format='png')
```

Plots the distribution of per gene minimum stabilities across all dataset pairs and returns a sorted dataframe of minimum gene stabilities. See `run_second_order_correlation_analysis` for more info.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `gene_stabilities` | np.ndarray | Stabilities tensor of shape `(n_datasets, n_datasets, n_genes)`. |
| `robust_gene_names` | np.ndarray | Gene names corresponding to the third dimension of `gene_stabilities`. |
| `node_label_df` | pd.DataFrame, optional | Node labels DataFrame with `'name'` and `'module_label'` columns. If provided, module labels are included in output. |
| `output_folder` | str, optional | Directory to save the plot. If `None`, does not save. |
| `figsize` | tuple, optional | Figure size in inches. Default is `(3, 2)`. |
| `fontsize` | int, optional | Font size for axes labels. Default is 6. |
| `bins` | int, optional | Number of histogram bins. Default is 50. |
| `dpi` | int, optional | Resolution for saved figure. Default is 300. |
| `file_format` | str, optional | `'png'` or `'pdf'`. Default is `'png'`. |

**Returns:** pd.DataFrame with columns `['gene', 'module', 'min_stability']` sorted from most to least stable. If `node_label_df` is not provided, `'module'` column will be NaN.

---

### `plot_gene_stability`

```python
plot_gene_stability(dataset_names, gene_name, gene_stabilities, robust_gene_names,
                    output_folder=None, figsize=(1.5, 1.5), fontsize=5,
                    dpi=300, file_format='png', x_ticks=True, cbar=True)
```

Plots a compact stability heatmap for a specific gene showing pairwise dataset correlations. See `run_second_order_correlation_analysis` for more info.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `dataset_names` | list of str | Dataset names for axis labels. |
| `gene_name` | str | Name of the gene to visualize. |
| `gene_stabilities` | np.ndarray | Stabilities tensor of shape `(n_datasets, n_datasets, n_genes)`. |
| `robust_gene_names` | np.ndarray | Gene names corresponding to the third dimension of `gene_stabilities`. |
| `output_folder` | str, optional | Directory to save the plot. If `None`, does not save. |
| `figsize` | tuple, optional | Figure size in inches. Default is `(1.5, 1.5)`. |
| `fontsize` | int, optional | Font size for title and tick labels. Default is 5. |
| `dpi` | int, optional | Resolution for saved figure. Default is 300. |
| `file_format` | str, optional | `'png'` or `'pdf'`. Default is `'png'`. |
| `x_ticks` | bool, optional | Whether to display x-axis dataset labels. Default is `True`. |
| `cbar` | bool, optional | Whether to display a colorbar. Default is `True`. |

**Returns:** None.

---

### `plot_module_stability`

```python
plot_module_stability(module_label, gene_stabilities, robust_gene_names, node_label_df,
                      dataset_names, output_folder=None, figsize=(1.5, 1.5), fontsize=5,
                      dpi=300, file_format='png', x_ticks=False, cbar=True)
```

Plots a compact stability heatmap for a module showing average pairwise dataset correlations across all genes in the module. See `run_second_order_correlation_analysis` for more info.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `module_label` | int or str | Module identifier from `node_label_df`. |
| `gene_stabilities` | np.ndarray | Stabilities tensor of shape `(n_datasets, n_datasets, n_genes)`. |
| `robust_gene_names` | np.ndarray | Gene names corresponding to the third dimension of `gene_stabilities`. |
| `node_label_df` | pd.DataFrame | DataFrame with `'name'` and `'module_label'` columns. |
| `dataset_names` | list of str | Dataset names for axis labels. |
| `output_folder` | str, optional | Directory to save the plot. If `None`, does not save. |
| `figsize` | tuple, optional | Figure size in inches. Default is `(1.5, 1.5)`. |
| `fontsize` | int, optional | Font size for title and labels. Default is 5. |
| `dpi` | int, optional | Resolution for saved figure. Default is 300. |
| `file_format` | str, optional | `'png'` or `'pdf'`. Default is `'png'`. |
| `x_ticks` | bool, optional | Whether to display x-axis dataset labels. Default is `False`. |
| `cbar` | bool, optional | Whether to display a colorbar. Default is `True`. |

**Returns:** None.

---

## choosing_hyperparameters.py

### `select_clustering_params`

```python
select_clustering_params(gene_names, pearsonR_mat, permuted_pcc_999=None, output_folder=None,
                          pcc_cutoffs=[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
                          clustering_powers=[1, 3, 5, 7, 9], min_genes_for_module=3,
                          infomap_clustering_trials=1, full_metrics=False)
```

Evaluates combinations of PCC cutoff and clustering power hyperparameters for spatial network construction, generating diagnostic plots of clustering quality metrics to guide hyperparameter selection.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `gene_names` | list, np.ndarray, or pd.Index | Gene names corresponding to rows/columns of `pearsonR_mat`. |
| `pearsonR_mat` | np.ndarray | Square matrix of Pearson correlation coefficients between genes. |
| `permuted_pcc_999` | float, optional | The 99.9th percentile PCC from a shuffled-data null distribution (output of `compute_shuffled_correlation_percentiles`). If provided, appended to `pcc_cutoffs` as an additional cutoff to evaluate. Default is `None`. |
| `output_folder` | str, optional | Directory where output plots and CSV results are saved. If `None`, plots are displayed but not saved. Default is `None`. |
| `pcc_cutoffs` | list of float, optional | PCC hard-threshold cutoff values to evaluate. Default is `[0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]`. |
| `clustering_powers` | list of int or float, optional | Soft-power exponents to evaluate for rescaling PCC edge weights. Default is `[1, 3, 5, 7, 9]`. |
| `min_genes_for_module` | int, optional | Minimum number of genes required for a cluster to be counted as a module. Default is 3. |
| `infomap_clustering_trials` | int, optional | Number of Infomap algorithm trials per hyperparameter combination. More trials improve stability. Default is 1 for runtime efficiency. |
| `full_metrics` | bool, optional | If `True`, plots all metrics (`mean_gene_margin`, `fraction_margin_positive`, `modularity`, `n_clusters`, `n_genes_included`). If `False`, plots only `mean_gene_margin`, `n_clusters`, and `n_genes_included`. Default is `False`. |

**Returns:** None. Displays diagnostic plots and optionally saves them along with a results CSV to `output_folder`.

---

## plotting.py

### `rotate_spatial`

```python
rotate_spatial(adata, angle_degrees=0, flip_x=False, flip_y=False,
               spatial_key='spatial', center=True)
```

Flips and rotates the spatial coordinates in an AnnData object.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `adata` | AnnData | The AnnData object to modify. |
| `angle_degrees` | float, optional | Angle to rotate in degrees (counter-clockwise). Default is 0. |
| `flip_x` | bool, optional | If `True`, mirrors the tissue horizontally. Default is `False`. |
| `flip_y` | bool, optional | If `True`, mirrors the tissue vertically. Default is `False`. |
| `spatial_key` | str, optional | Key in `adata.obsm` storing `(x, y)` coordinates. Default is `'spatial'`. |
| `center` | bool, optional | If `True`, flips and rotates around the tissue center. If `False`, rotates around the `(0,0)` origin. Default is `True`. |

**Returns:** AnnData — The updated AnnData object.

---

### `plot_gene`

```python
plot_gene(sm_adata, gene_name, output_folder=None, figsize=(1,1), spot_size=25,
          fontsize=6, fontfamily='sans-serif', cmap=None, dpi=300, file_format='png')
```

Plots the spatial expression of a specified gene in a single spatial transcriptomic dataset.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `sm_adata` | AnnData | AnnData object containing spatial transcriptomic data. |
| `gene_name` | str | Name of the gene to visualize. |
| `output_folder` | str, optional | Directory path where the output plot is saved. If `None`, displayed but not saved. Default is `None`. |
| `figsize` | tuple, optional | Figure size in inches `(width, height)`. Default is `(1, 1)`. |
| `spot_size` | int, optional | Size of spatial spots in the plot. Default is 25. |
| `fontsize` | int, optional | Font size for the plot title. Default is 6. |
| `fontfamily` | str, optional | Font family for the plot title. Default is `'sans-serif'`. |
| `cmap` | Colormap, optional | Colormap for visualizing expression. If `None`, defaults to a custom gray-red-black colormap. |
| `dpi` | int, optional | Resolution of the saved plot in DPI. Default is 300. |
| `file_format` | str, optional | `'png'` or `'pdf'`. Default is `'png'`. |

**Returns:** None.

---

### `plot_modules`

```python
plot_modules(sm_adata, node_label_df, output_folder=None, plots_per_row=5, min_genes=3,
             figsize=None, spot_size=25, fontsize=6, fontfamily='sans-serif',
             cmap=None, dpi=300, file_format='png')
```

Plots module scores spatially, arranging plots in rows with a fixed number of columns.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `sm_adata` | AnnData | AnnData object containing spatial transcriptomic data. |
| `node_label_df` | pd.DataFrame | Gene module assignments with columns `'module_label'` and `'name'`. |
| `output_folder` | str, optional | Directory path where output plots are saved. If `None`, displayed but not saved. Default is `None`. |
| `plots_per_row` | int, optional | Number of plots per row in each output file. Default is 5. |
| `min_genes` | int, optional | Minimum number of genes required for a module to be plotted. Default is 3. |
| `figsize` | tuple, optional | Figure size in inches. Default is `(plots_per_row, 1)`. |
| `spot_size` | int, optional | Size of spatial spots. Default is 25. |
| `fontsize` | int, optional | Font size for subplot titles. Default is 6. |
| `fontfamily` | str, optional | Font family for plot titles. Default is `'sans-serif'`. |
| `cmap` | Colormap, optional | Colormap for visualizing expression. If `None`, defaults to a custom gray-red-black colormap. |
| `dpi` | int, optional | Resolution of saved plots in DPI. Default is 300. |
| `file_format` | str, optional | `'png'` or `'pdf'`. Default is `'png'`. |

**Returns:** None.

---

### `plot_gene_multisample`

```python
plot_gene_multisample(sm_adata_list, adata_list_names, gene_name, output_folder=None,
                      shared_scaling=False, figsize=None, spot_size=25, fontsize=6,
                      fontfamily='sans-serif', cmap=None, dpi=300, file_format='png')
```

Plots the spatial expression of a specified gene across multiple spatial transcriptomic datasets.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `sm_adata_list` | list of AnnData | AnnData objects for each dataset. |
| `adata_list_names` | list of str | Dataset names for labeling subplots. |
| `gene_name` | str | Name of the gene to visualize. |
| `output_folder` | str, optional | Directory where output plots are saved. If `None`, displayed but not saved. Default is `None`. |
| `shared_scaling` | bool, optional | If `True`, normalizes using the maximum expression value across all datasets. If `False`, normalizes within each dataset independently. Default is `False`. |
| `figsize` | tuple, optional | Figure size in inches. Default is `(n_datasets, 1)`. |
| `spot_size` | int, optional | Size of spatial spots. Default is 25. |
| `fontsize` | int, optional | Font size for subplot titles. Default is 6. |
| `fontfamily` | str, optional | Font family for text annotations. Default is `'sans-serif'`. |
| `cmap` | Colormap, optional | Colormap for visualizing expression. If `None`, defaults to a custom gray-red-black colormap. |
| `dpi` | int, optional | Resolution of saved plots in DPI. Default is 300. |
| `file_format` | str, optional | `'png'` or `'pdf'`. Default is `'png'`. |

**Returns:** None.

---

### `plot_modules_multisample`

```python
plot_modules_multisample(sm_adata_list, adata_list_names, node_label_df, output_folder=None,
                          shared_scaling=False, min_genes=3, figsize=None, spot_size=25,
                          fontsize=6, fontfamily='sans-serif', cmap=None, dpi=300, file_format='png')
```

Plots gene module scores across multiple spatial transcriptomic datasets, with an option for shared scaling.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `sm_adata_list` | list of AnnData | AnnData objects for each dataset. |
| `adata_list_names` | list of str | Dataset names for labeling subplots. |
| `node_label_df` | pd.DataFrame | Gene module assignments with columns `'module_label'` and `'name'`. |
| `output_folder` | str, optional | Directory where output plots are saved. If `None`, displayed but not saved. Default is `None`. |
| `shared_scaling` | bool, optional | If `True`, normalizes using the maximum expression value across all datasets. If `False`, normalizes within each dataset independently. Default is `False`. |
| `min_genes` | int, optional | Minimum number of genes required for a module to be plotted. Default is 3. |
| `figsize` | tuple, optional | Figure size in inches. Default is `(n_datasets, 1)`. |
| `spot_size` | int, optional | Size of spatial spots. Default is 25. |
| `fontsize` | int, optional | Font size for subplot titles. Default is 6. |
| `fontfamily` | str, optional | Font family for text annotations. Default is `'sans-serif'`. |
| `cmap` | Colormap, optional | Colormap for visualizing expression. If `None`, defaults to a custom gray-red-black colormap. |
| `dpi` | int, optional | Resolution of saved plots in DPI. Default is 300. |
| `file_format` | str, optional | `'png'` or `'pdf'`. Default is `'png'`. |

**Returns:** None.

---

## create_anndata_from_transcripts.py

### `create_anndata_from_transcripts`

```python
create_anndata_from_transcripts(input_file, x_col, y_col, gene_col, output_file,
                                 file_format='csv', min_counts_per_gene=1,
                                 count_col=None, chunksize=1000000)
```

Converts submicron-level spatial data (MERFISH, Xenium, CosMx, seqFISH, Stereo-seq bin1) to AnnData format. Each unique coordinate becomes a separate spot, preserving coordinates in their original units.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `input_file` | str | Path to input file (CSV, Parquet, TSV, or TSV.GZ). |
| `x_col` | str | Column name for x coordinates (e.g., `'x_location'`, `'global_x'`, `'x'`). |
| `y_col` | str | Column name for y coordinates (e.g., `'y_location'`, `'global_y'`, `'y'`). |
| `gene_col` | str | Column name for gene names (e.g., `'feature_name'`, `'gene'`, `'target'`, `'geneID'`). |
| `output_file` | str | Path to save output `.h5ad` file. If `None`, does not save to disk. |
| `file_format` | str, optional | Input file format: `'csv'`, `'parquet'`, `'tsv'`, or `'tsv.gz'`. Default is `'csv'`. |
| `min_counts_per_gene` | int, optional | Minimum number of transcripts per gene to include in output. Default is 1. |
| `count_col` | str, optional | Column name for pre-aggregated counts (e.g., `'MIDCounts'`). If `None`, assumes one transcript per row. |
| `chunksize` | int, optional | Number of rows to read at a time for TSV/TSV.GZ files. Default is 1,000,000. |

**Returns:** AnnData with:

- `X`: Sparse CSR matrix of gene counts (spots × genes)
- `obs`: Spot metadata with `'total_counts'`
- `var`: Gene metadata with `'total_counts'`
- `obsm['spatial']`: Spatial coordinates `(n_spots × 2)`

**Example Usage:**

> **Command line recommended** — run in a tmux screen or equivalent for long-running conversions.

**Command Line:**

Xenium (10x Genomics):
```bash
python -m smoothie.create_anndata_from_transcripts \
    transcripts.parquet \
    --x-col x_location \
    --y-col y_location \
    --gene-col feature_name \
    --format parquet \
    --output xenium_data.h5ad
```

Stereo-seq (MGI/BGI):
```bash
python -m smoothie.create_anndata_from_transcripts \
    E9.5_E1S1_GEM_bin1.tsv.gz \
    --x-col x \
    --y-col y \
    --gene-col geneID \
    --count-col MIDCounts \
    --format tsv.gz \
    --chunksize 1000000 \
    --output stereoseq_data.h5ad
```

MERFISH (Vizgen):
```bash
python -m smoothie.create_anndata_from_transcripts \
    detected_transcripts.csv \
    --x-col global_x \
    --y-col global_y \
    --gene-col gene \
    --format csv \
    --output merfish_data.h5ad
```

CosMx (NanoString):
```bash
python -m smoothie.create_anndata_from_transcripts \
    transcripts.csv \
    --x-col x_global_px \
    --y-col y_global_px \
    --gene-col target \
    --format csv \
    --output cosmx_data.h5ad
```

seqFISH/seqFISH+:
```bash
python -m smoothie.create_anndata_from_transcripts \
    transcripts.csv \
    --x-col x \
    --y-col y \
    --gene-col gene \
    --format csv \
    --output seqfish_data.h5ad
```

**Python:**

Xenium (10x Genomics):
```python
adata = smoothie.create_anndata_from_transcripts(
    'transcripts.parquet',
    x_col='x_location',
    y_col='y_location',
    gene_col='feature_name',
    file_format='parquet',
    output_file='xenium_data.h5ad'
)
```

Stereo-seq (MGI/BGI):
```python
adata = smoothie.create_anndata_from_transcripts(
    'E9.5_E1S1_GEM_bin1.tsv.gz',
    x_col='x',
    y_col='y',
    gene_col='geneID',
    count_col='MIDCounts',
    file_format='tsv.gz',
    output_file='stereoseq_data.h5ad'
)
```

MERFISH (Vizgen):
```python
adata = smoothie.create_anndata_from_transcripts(
    'detected_transcripts.csv',
    x_col='global_x',
    y_col='global_y',
    gene_col='gene',
    file_format='csv',
    output_file='merfish_data.h5ad'
)
```

CosMx (NanoString):
```python
adata = smoothie.create_anndata_from_transcripts(
    'transcripts.csv',
    x_col='x_global_px',
    y_col='y_global_px',
    gene_col='target',
    file_format='csv',
    output_file='cosmx_data.h5ad'
)
```

seqFISH/seqFISH+:
```python
adata = smoothie.create_anndata_from_transcripts(
    'transcripts.csv',
    x_col='x',
    y_col='y',
    gene_col='gene',
    file_format='csv',
    output_file='seqfish_data.h5ad'
)
```

---

## shuffle_analysis.py

### `compute_shuffled_correlation_percentiles`

```python
compute_shuffled_correlation_percentiles(adata, grid_based_or_not, gaussian_sd=20,
                                          stride=None, grid_fitting_dist=None,
                                          min_spots_under_gaussian=25, num_processes=4,
                                          num_data_splits=None, seed=None)
```

Shuffles spatial coordinates, smoothes data, and computes correlation matrix percentiles for null distribution estimation.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `adata` | AnnData or list of AnnData | Single AnnData object or list of AnnData objects. For multiple datasets, genes are aligned. |
| `grid_based_or_not` | bool | `True` = bin shuffling + grid-based smoothing (for subcellular resolution data); `False` = in-place shuffling + in-place smoothing (for cellular resolution data). |
| `gaussian_sd` | float | Standard deviation of Gaussian kernel in coordinate units. |
| `min_spots_under_gaussian` | int | Minimum spots required within `3 * gaussian_sd` radius for valid smoothing at a given location. |
| `stride` | float, optional | Grid spacing for smoothing (grid-based only). Default is `gaussian_sd`. |
| `grid_fitting_dist` | float, optional | Minimum distance from grid point to data for grid fitting (grid-based only). Default is `0.25 × gaussian_sd`. |
| `num_processes` | int, optional | Number of parallel processes for smoothing. Default is 4. |
| `num_data_splits` | int, optional | Number of data splits for parallel processing. If `None`, automatically determined. |
| `seed` | int, optional | Random seed for reproducibility. Default is `None`. |

**Returns:** `(p95, p99, p999)` — The 95th, 99th, and 99.9th percentiles of the shuffled correlation coefficient distribution. Use these as `pcc_cutoff` inputs to `select_clustering_params` or `make_spatial_network`.

> **Note:** Call this function with the same parameters used for `run_parallelized_smoothing`.

---

## utils.py

### `suppress_warnings`

```python
suppress_warnings()
```

Suppresses specific warnings that Smoothie commonly triggers due to memory-efficient operations, including AnnData `ImplicitModificationWarning` and scanpy view warnings. Called automatically on import.

**Returns:** None.

---

### `enable_warnings`

```python
enable_warnings()
```

Re-enables all warnings that were suppressed by resetting all warning filters to default.

**Returns:** None.

---

### `quiet_mode`

```python
quiet_mode()
```

Context manager for temporarily suppressing Smoothie warnings in a specific code block. Warnings are automatically restored after the block exits.

**Returns:** Context manager.

**Example:**
```python
with smoothie.quiet_mode():
    sm_adata = smoothie.run_parallelized_smoothing(...)
# warnings restored after the block
```

<img src="../images/hidden_img.png" alt="surprise" width="500"/>