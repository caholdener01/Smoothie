Smoothie expects your spatial count matrix to be in the AnnData `.h5ad` format. Here are quick instructions on creating the AnnData object from different spatial transcriptomics platforms outputs. See AnnData documentation here: [AnnData Docs](https://anndata.readthedocs.io/en/stable/).

### **Option 1: Use Scanpy to load count matrices**
#### Slide-seq
```python
import scanpy as sc
import pandas as pd

# Update these paths accordingly
adata = sc.read_10x_mtx("path_to_mtx/")
coords = pd.read_csv("bead_locations.csv")
adata.obsm["spatial"] = coords[["x", "y"]].values

# save output
adata.write_h5ad("dataset.h5ad")
```
#### VisiumHD
```python
import numpy as np
import pandas as pd

# Update these paths accordingly
adata = sc.read_10x_mtx("./binned_outputs/square_002um/filtered_feature_bc_matrix/")
coords = pd.read_parquet("./binned_outputs/square_002um/spatial/tissue_positions.parquet")

# Ensure barcode column is string (important!)
coords['barcode'] = coords['barcode'].astype(str)
adata.obs_names = adata.obs_names.astype(str)

# Set barcode as index and reindex to match AnnData order
coords_indexed = coords.set_index('barcode')
coords_matched = coords_indexed.loc[adata.obs_names]
assert np.all(coords_matched.index == adata.obs_names)

# Copy spatial coordinates
adata.obsm["spatial"] = coords_matched[
    ["pxl_col_in_fullres", "pxl_row_in_fullres"]
].values

# save output
adata.write_h5ad("dataset.h5ad")
```

### **Option 2: Use smoothie.create_anndata_from_transcripts**
This works for most imaging-based platforms and Stereo-seq bin1 output. Runtime can take minutes to hours, so *I recommend running this code in the terminal in the background with tmux screens or an equivalent*.

#### Xenium, MERFISH, CosMx, seqFISH, Stereo-Seq (bin1)

Python
```python
import smoothie
adata = smoothie.create_anndata_from_transcripts(
    input_file='transcripts.parquet', # path to transcripts file 
    x_col='x_location', # name of x coordinate column
    y_col='y_location', # name of y coordinate column
    gene_col='feature_name', # name of gene column
    file_format='parquet', # 'csv', 'parquet', 'tsv', or 'tsv.gz'
    output_file='dataset.h5ad' # path to save AnnData .h5ad
)
```
Bash (Terminal)
```bash
python -m smoothie.create_anndata_from_transcripts \
        transcripts.parquet \
        --x-col x_location \
        --y-col y_location \
        --gene-col feature_name \
        --format parquet \
        --output xenium_data.h5ad
```
Note, the code above is configured for Xenium data. Check the API for **smoothie.create_anndata_from_transcripts** for settings to use for MERFISH, CosMx, seqFish, Stereo-seq (bin1).

### **Option 3: Build from Components**
This works if you are moving from R to Python.
```python
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix

# load count_matrix, gene_names_arr, coords with pandas

# ensure count matrix is CSR sparse format
count_matrix = csr_matrix(count_matrix) # spots (rows) x genes (cols)

adata = ad.AnnData(X=count_matrix)
adata.var_names=gene_names_arr # gene names for each column of count matrix
adata.obsm["spatial"] = coords[["x", "y"]].values

# save output
adata.write_h5ad("dataset.h5ad")
```