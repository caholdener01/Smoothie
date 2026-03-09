# Copyright (c) 2026 Chase Holdener
# Licensed under the MIT License. See LICENSE file for details.

"""
Utility for converting imaging-based spatial transcriptomics data to AnnData format.

This module provides functions to convert submicron resolution CSV or Parquet files from 
imaging-based spatial transcriptomics platforms (e.g., MERFISH, seqFISH, Xenium, CosMx)
or Stereo-seq bin1.tsv.gz files into AnnData objects compatible with Smoothie analysis.
"""

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix, dok_matrix
from pathlib import Path
import pyarrow.parquet as pq
import warnings


def create_anndata_from_transcripts(
    input_file: str,
    x_col: str,
    y_col: str,
    gene_col: str,
    output_file: str,
    file_format: str = 'csv',
    min_counts_per_gene: int = 1,
    count_col: str = None,
    chunksize: int = 1000000
) -> ad.AnnData:
    """
    Convert "transcript-level" spatial data (MERFISH, Xenium, CosMx, seqFISH, bin1 Stereo-seq) to AnnData format.
    
    This function processes imaging-based spatial transcriptomics data at native resolution,
    where each unique coordinate becomes a separate spot. Coordinates are preserved in their
    original units (microns, pixels, etc.) with no conversion applied.
    
    Parameters
    ----------
    input_file : str
        Path to input file (CSV, Parquet, TSV, or TSV.GZ format).
    x_col : str
        Column name for x coordinates (e.g., 'x_location', 'global_x', 'x').
    y_col : str
        Column name for y coordinates (e.g., 'y_location', 'global_y', 'y').
    gene_col : str
        Column name for gene names (e.g., 'feature_name', 'gene', 'target', 'geneID').
    output_file : str, optional
        Path to save output .h5ad file. If None, does not save to disk.
    file_format : str, default='csv'
        Format of input file: 'csv', 'parquet', 'tsv', or 'tsv.gz'.
    min_counts_per_gene : int, default=1
        Minimum number of transcripts per gene to include in output.
    count_col : str, optional
        Column name for pre-aggregated counts (e.g., 'MIDCounts', 'count').
        If None, assumes one transcript per row.
    chunksize : int, default=1000000
        Number of rows to read at a time for TSV/TSV.GZ files (memory efficiency).
    
    Returns
    -------
    ad.AnnData
        AnnData object with:
        - X: Sparse CSR matrix of gene counts (spots × genes)
        - obs: Spot metadata with 'total_counts'
        - var: Gene metadata with 'total_counts'
        - obsm['spatial']: Spatial coordinates (n_spots × 2)
    
    Examples
    --------
    COMMAND LINE: Platform-specific examples 
    (** Recommended: Run command in the background (with tmux screen or equivalent):
    
    Xenium (10x Genomics):
    $ python -m smoothie.create_anndata_from_transcripts \
        transcripts.parquet \
        --x-col x_location \
        --y-col y_location \
        --gene-col feature_name \
        --format parquet \
        --output xenium_data.h5ad
    
    Stereo-seq (MGI/BGI):
    $ python -m smoothie.create_anndata_from_transcripts \
        E9.5_E1S1_GEM_bin1.tsv.gz \
        --x-col x \
        --y-col y \
        --gene-col geneID \
        --count-col MIDCounts \
        --format tsv.gz \
        --chunksize 1000000 \
        --output stereoseq_data.h5ad
    
    MERFISH (Vizgen):
    $ python -m smoothie.create_anndata_from_transcripts \
        detected_transcripts.csv \
        --x-col global_x \
        --y-col global_y \
        --gene-col gene \
        --format csv \
        --output merfish_data.h5ad
    
    CosMx (NanoString):
    $ python -m smoothie.create_anndata_from_transcripts \
        transcripts.csv \
        --x-col x_global_px \
        --y-col y_global_px \
        --gene-col target \
        --format csv \
        --output cosmx_data.h5ad
    
    seqFISH/seqFISH+:
    $ python -m smoothie.create_anndata_from_transcripts \
        transcripts.csv \
        --x-col x \
        --y-col y \
        --gene-col gene \
        --format csv \
        --output seqfish_data.h5ad
    
    
    PYTHON: Platform-specific examples
    (** Note: Long runtimes... I recommend running in command line instead (See above usage):
    
    Xenium (10x Genomics) - Parquet format:
    >>> adata = create_anndata_from_transcripts(
    ...     'transcripts.parquet',
    ...     x_col='x_location',
    ...     y_col='y_location', 
    ...     gene_col='feature_name',
    ...     file_format='parquet',
    ...     output_file='xenium_data.h5ad'
    ... )
    
    Stereo-seq - TSV.GZ with pre-aggregated counts:
    >>> adata = create_anndata_from_transcripts(
    ...     'E9.5_E1S1_GEM_bin1.tsv.gz',
    ...     x_col='x',
    ...     y_col='y',
    ...     gene_col='geneID',
    ...     count_col='MIDCounts',
    ...     file_format='tsv.gz',
    ...     output_file='stereoseq_data.h5ad'
    ... )
    
    MERFISH (Vizgen) - CSV format:
    >>> adata = create_anndata_from_transcripts(
    ...     'detected_transcripts.csv',
    ...     x_col='global_x',
    ...     y_col='global_y',
    ...     gene_col='gene',
    ...     file_format='csv',
    ...     output_file='merfish_data.h5ad'
    ... )
    
    CosMx (NanoString) - CSV format:
    >>> adata = create_anndata_from_transcripts(
    ...     'transcripts.csv',
    ...     x_col='x_global_px',
    ...     y_col='y_global_px',
    ...     gene_col='target',
    ...     file_format='csv',
    ...     output_file='cosmx_data.h5ad'
    ... )
    
    seqFISH/seqFISH+ - CSV format:
    >>> adata = create_anndata_from_transcripts(
    ...     'transcripts.csv',
    ...     x_col='x',
    ...     y_col='y',
    ...     gene_col='gene',
    ...     file_format='csv',
    ...     output_file='seqfish_data.h5ad'
    ... )
    
    Notes
    -----
    - Native resolution: Each unique (x,y) coordinate becomes a separate spot
    - Coordinates are preserved in their original scale (no unit conversion)
    - Floating-point coordinates are rounded to 6 decimal places for grouping
    - Output matrix is automatically sparse (CSR format) for memory efficiency
    - Memory usage is reported after conversion
    - Expect high sparsity (>99%) which is normal and efficient
    - TSV.GZ files are read in chunks for memory efficiency with large datasets
    - If you need binning, use scanpy's spatial binning functions after conversion
    """

    if output_file:
        out_dir = Path(output_file).parent
        # If out_dir isn't just the current directory (empty string path) and doesn't exist
        if out_dir.name and not out_dir.exists():
            raise FileNotFoundError(
                f"Output directory does not exist: '{out_dir}'. "
                f"Please create it before running the script."
            )
    
    print(f"Reading {file_format.upper()} file: {input_file}")
    
    # Determine file reading parameters
    if file_format.lower() == 'tsv.gz':
        sep = '\t'
        compression = 'gzip'
        use_chunks = True
    elif file_format.lower() == 'tsv':
        sep = '\t'
        compression = None
        use_chunks = True
    elif file_format.lower() == 'csv':
        sep = ','
        compression = None
        use_chunks = True
    elif file_format.lower() == 'parquet':
        # Parquet requires different handling - read in batches
        use_chunks = False
    else:
        raise ValueError(f"Unsupported file format: {file_format}. Use 'csv', 'parquet', 'tsv', or 'tsv.gz'.")
    
    if use_chunks:
        # CHUNKED READING PATH (CSV, TSV, TSV.GZ)
        print(f"  Using chunked reading (chunksize={chunksize:,} rows)")
        
        # Step 1: First pass - collect unique coordinates and genes
        print("  Pass 1/2: Collecting unique coordinates and genes...")
        xy_set = set()
        gene_set = set()
        
        reader = pd.read_csv(input_file, sep=sep, compression=compression, chunksize=chunksize)
        
        chunk_count = 0
        total_rows = 0
        for chunk in reader:
            # Validate columns on first chunk
            if chunk_count == 0:
                required_cols = [x_col, y_col, gene_col]
                if count_col:
                    required_cols.append(count_col)
                missing = [col for col in required_cols if col not in chunk.columns]
                if missing:
                    raise ValueError(f"Missing columns in input file: {missing}\nAvailable columns: {chunk.columns.tolist()}")
            
            xy_set.update(set(zip(chunk[x_col].round(6), chunk[y_col].round(6))))
            gene_set.update(chunk[gene_col].dropna().unique())
            total_rows += len(chunk)
            chunk_count += 1
            if chunk_count % 10 == 0:
                print(f"    Processed {total_rows:,} rows...")
        
        print(f"  Total rows: {total_rows:,}")
        print(f"  Found {len(xy_set):,} unique coordinates, {len(gene_set):,} genes")
        
        # Create mappings
        xy_list = sorted(list(xy_set))
        gene_list = sorted(list(gene_set))
        xy_dict = {xy: i for i, xy in enumerate(xy_list)}
        gene_dict = {gene: i for i, gene in enumerate(gene_list)}
        
        # Filter genes by minimum counts if needed
        if min_counts_per_gene > 1:
            print(f"  Counting gene totals for filtering (min_counts={min_counts_per_gene})...")
            gene_totals = {gene: 0 for gene in gene_list}
            reader = pd.read_csv(input_file, sep=sep, compression=compression, chunksize=chunksize)
            for chunk in reader:
                if count_col:
                    gene_chunk_counts = chunk.groupby(gene_col)[count_col].sum()
                else:
                    gene_chunk_counts = chunk[gene_col].value_counts()
                for gene, count in gene_chunk_counts.items():
                    if gene in gene_totals:
                        gene_totals[gene] += count
            
            valid_genes = [gene for gene, count in gene_totals.items() if count >= min_counts_per_gene]
            print(f"  Filtered genes: {len(valid_genes)}/{len(gene_list)} genes kept (≥{min_counts_per_gene} counts)")
            
            gene_list = sorted(valid_genes)
            gene_dict = {gene: i for i, gene in enumerate(gene_list)}
        
        # Initialize sparse matrix
        num_rows = len(xy_list)
        num_cols = len(gene_list)
        print(f"  Initializing sparse matrix: {num_rows:,} spots × {num_cols:,} genes")
        sparse_matrix = dok_matrix((num_rows, num_cols), dtype=np.float32)
        
        # Step 2: Second pass - fill sparse matrix
        print("  Pass 2/2: Building count matrix...")
        reader = pd.read_csv(input_file, sep=sep, compression=compression, chunksize=chunksize)
        
        chunk_count = 0
        rows_processed = 0
        for chunk in reader:
            # Filter to valid genes
            chunk = chunk[chunk[gene_col].isin(gene_dict.keys())]
            
            for idx, row in chunk.iterrows():
                x = row[x_col]
                y = row[y_col]
                gene = row[gene_col]
                count = row[count_col] if count_col else 1
                
                xy_idx = xy_dict.get((x, y))
                gene_idx = gene_dict.get(gene)
                
                if xy_idx is not None and gene_idx is not None:
                    sparse_matrix[xy_idx, gene_idx] += count
            
            rows_processed += len(chunk)
            chunk_count += 1
            if chunk_count % 10 == 0:
                print(f"    Processed {rows_processed:,} rows...")
        
        # Convert to CSR format
        print("  Converting to CSR format...")
        sparse_matrix = sparse_matrix.tocsr()
        
        # Get coordinates
        spatial_coords = np.array(xy_list, dtype=np.float32)
        
        print(f"  Total counts in matrix: {int(sparse_matrix.sum()):,}")
    
    else:
        # PARQUET PATH (load in batches using pyarrow)
        print(f"  Reading Parquet file in batches...")
        
        parquet_file = pq.ParquetFile(input_file)
        total_rows = parquet_file.metadata.num_rows
        print(f"  Total rows: {total_rows:,}")
        
        # Validate columns
        available_cols = parquet_file.schema.names
        required_cols = [x_col, y_col, gene_col]
        if count_col:
            required_cols.append(count_col)
        missing = [col for col in required_cols if col not in available_cols]
        if missing:
            raise ValueError(f"Missing columns in Parquet file: {missing}\nAvailable columns: {available_cols}")
        
        # Step 1: First pass - collect unique coordinates and genes
        print("  Pass 1/2: Collecting unique coordinates and genes...")
        xy_set = set()
        gene_set = set()
        
        batch_count = 0
        for batch in parquet_file.iter_batches(batch_size=chunksize, columns=[x_col, y_col, gene_col]):
            df_batch = batch.to_pandas()
            xy_set.update(set(zip(df_batch[x_col], df_batch[y_col])))
            gene_set.update(df_batch[gene_col].dropna().unique())
            batch_count += 1
            if batch_count % 10 == 0:
                print(f"    Processed {batch_count * chunksize:,} rows...")
        
        print(f"  Found {len(xy_set):,} unique coordinates, {len(gene_set):,} genes")
        
        # Create mappings
        xy_list = sorted(list(xy_set))
        gene_list = sorted(list(gene_set))
        xy_dict = {xy: i for i, xy in enumerate(xy_list)}
        gene_dict = {gene: i for i, gene in enumerate(gene_list)}
        
        # Filter genes by minimum counts if needed
        if min_counts_per_gene > 1:
            print(f"  Counting gene totals for filtering (min_counts={min_counts_per_gene})...")
            gene_totals = {gene: 0 for gene in gene_list}
            
            cols_to_read = [gene_col, count_col] if count_col else [gene_col]
            for batch in parquet_file.iter_batches(batch_size=chunksize, columns=cols_to_read):
                df_batch = batch.to_pandas()
                if count_col:
                    gene_chunk_counts = df_batch.groupby(gene_col)[count_col].sum()
                else:
                    gene_chunk_counts = df_batch[gene_col].value_counts()
                for gene, count in gene_chunk_counts.items():
                    if gene in gene_totals:
                        gene_totals[gene] += count
            
            valid_genes = [gene for gene, count in gene_totals.items() if count >= min_counts_per_gene]
            print(f"  Filtered genes: {len(valid_genes)}/{len(gene_list)} genes kept (≥{min_counts_per_gene} counts)")
            
            gene_list = sorted(valid_genes)
            gene_dict = {gene: i for i, gene in enumerate(gene_list)}
        
        # Initialize sparse matrix
        num_rows = len(xy_list)
        num_cols = len(gene_list)
        print(f"  Initializing sparse matrix: {num_rows:,} spots × {num_cols:,} genes")
        sparse_matrix = dok_matrix((num_rows, num_cols), dtype=np.float32)
        
        # Step 2: Second pass - fill sparse matrix
        print("  Pass 2/2: Building count matrix...")
        
        cols_to_read = [x_col, y_col, gene_col]
        if count_col:
            cols_to_read.append(count_col)
        
        batch_count = 0
        rows_processed = 0
        for batch in parquet_file.iter_batches(batch_size=chunksize, columns=cols_to_read):
            df_batch = batch.to_pandas()
            
            # Filter to valid genes
            df_batch = df_batch[df_batch[gene_col].isin(gene_dict.keys())]
            
            # Round coordinates
            df_batch[x_col] = df_batch[x_col].round(6)
            df_batch[y_col] = df_batch[y_col].round(6)
            
            for idx, row in df_batch.iterrows():
                x = row[x_col]
                y = row[y_col]
                gene = row[gene_col]
                count = row[count_col] if count_col else 1
                
                xy_idx = xy_dict.get((x, y))
                gene_idx = gene_dict.get(gene)
                
                if xy_idx is not None and gene_idx is not None:
                    sparse_matrix[xy_idx, gene_idx] += count
            
            rows_processed += len(df_batch)
            batch_count += 1
            if batch_count % 10 == 0:
                print(f"    Processed {rows_processed:,} rows...")
        
        # Convert to CSR format
        print("  Converting to CSR format...")
        sparse_matrix = sparse_matrix.tocsr()
        
        # Get coordinates
        spatial_coords = np.array(xy_list, dtype=np.float32)
        
        print(f"  Total counts in matrix: {int(sparse_matrix.sum()):,}")
    
    # Create AnnData object (common to both paths)
    print("  Creating AnnData object...")
    adata = ad.AnnData(X=sparse_matrix)
    adata.obsm['spatial'] = spatial_coords
    adata.var_names = gene_list
    
    # Save if requested
    if output_file:
        print(f"\nSaving to: {output_file}")
        adata.write_h5ad(output_file)
    
    return adata


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Convert transcript-level spatial data to AnnData format.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('input_file', type=str, help='Path to input file (CSV, Parquet, TSV, or TSV.GZ)')
    parser.add_argument('--x-col', required=True, type=str, help='Column name for x coordinates')
    parser.add_argument('--y-col', required=True, type=str, help='Column name for y coordinates')
    parser.add_argument('--gene-col', required=True, type=str, help='Column name for gene names')
    
    # Optional arguments
    parser.add_argument('--output', type=str, default=None, help='Path to save output .h5ad file')
    parser.add_argument('--format', type=str, default='csv', choices=['csv', 'parquet', 'tsv', 'tsv.gz'], help='Format of input file')
    parser.add_argument('--min-counts', type=int, default=1, help='Minimum number of transcripts per gene')
    parser.add_argument('--count-col', type=str, default=None, help='Column name for pre-aggregated counts')
    parser.add_argument('--chunksize', type=int, default=1000000, help='Number of rows to read at a time')

    args = parser.parse_args()

    # Run the function with the parsed arguments
    create_anndata_from_transcripts(
        input_file=args.input_file,
        x_col=args.x_col,
        y_col=args.y_col,
        gene_col=args.gene_col,
        output_file=args.output,
        file_format=args.format,
        min_counts_per_gene=args.min_counts,
        count_col=args.count_col,
        chunksize=args.chunksize
    )
