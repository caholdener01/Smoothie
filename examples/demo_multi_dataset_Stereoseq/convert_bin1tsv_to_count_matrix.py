# Chase Holdener, June 2024

print("Importing pandas...", end = '')
import pandas as pd
print("Importing numpy...", end = '')
import numpy as np
print("Importing anndata...", end = '')
import anndata as ad
from scipy.sparse import dok_matrix, csr_matrix


## CHOOSE PARAMETERS
FILE_NAME = './E9.5_E1S1_GEM_bin1.tsv.gz'
SAVE_SUFFIX_NAME = 'E9.5_E1S1'
CHUNKSIZE = 1000000


## Step 1: Read the file in chunks to find unique (x, y) pairs and gene values
xy_set = set()
gene_set = set()

reader = pd.read_csv(FILE_NAME, sep='\t', chunksize=CHUNKSIZE)

for chunk in reader:
    xy_set.update(set(zip(chunk['x'], chunk['y'])))
    gene_set.update(chunk['geneID'].unique())

xy_list = list(xy_set)
gene_list = list(gene_set)

xy_dict = {xy: i for i, xy in enumerate(xy_list)}
gene_dict = {gene: i for i, gene in enumerate(gene_list)}

# Initialize a sparse matrix with dok_matrix (can be converted to csr_matrix later)
num_rows = len(xy_list)
num_cols = len(gene_list)
sparse_matrix = dok_matrix((num_rows, num_cols), dtype=np.float32)


## Step 2: Fill the sparse matrix
reader = pd.read_csv(FILE_NAME, sep='\t', chunksize=CHUNKSIZE)

for chunk in reader:
    for idx, row in chunk.iterrows():
        x = row['x']
        y = row['y']
        gene = row['geneID']
        count = row['MIDCounts']
        
        # Get the row and column indices
        xy_idx = xy_dict[(x, y)]
        gene_idx = gene_dict[gene]
        
        # Update the sparse matrix
        sparse_matrix[xy_idx, gene_idx] += count

# Convert to CSR format for efficient arithmetic and matrix-vector operations
sparse_matrix = sparse_matrix.tocsr()


## Step 3: Get coordinates for each row of the sparse matrix

# Create reverse mappings for (x, y) indices
reverse_xy_dict = {v: k for k, v in xy_dict.items()}

# Create a DataFrame to store the row labels
row_labels = []
for row_idx in range(num_rows):
    x_value, y_value = reverse_xy_dict[row_idx]
    row_labels.append((x_value, y_value))

row_labels_df = pd.DataFrame(row_labels, columns=['x', 'y'])


## Step 4: Make adata and save it
dnb_coords = np.array(list(xy_dict.keys()))

embryo_adata = ad.AnnData(X = sparse_matrix)
embryo_adata.obs_names = [f"dnb_{i:d}" for i in range(embryo_adata.n_obs)]
embryo_adata.var_names = list(gene_dict.keys())
embryo_adata.obsm['spatial'] = dnb_coords

# Save anndata
embryo_adata.write_h5ad("mouse_"+SAVE_SUFFIX_NAME+"_raw_dnb.h5ad")
