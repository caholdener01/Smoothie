# Smoothie
Smoothie is a method that denoises spatial transcriptomics data with Gaussian smoothing and constructs and integrates genome-wide co-expression networks. Designed as a tool for biological discovery, Smoothie's output gene network allows for precise gene module detection, spatial annotation of uncharacterized genes, linkage of gene expression to genome architecture, and more. Further, Smoothie supports multi-sample comparisons to assess stable or dynamic gene expression patterns across tissues, conditions, and time points.  

## Cite Smoothie
"Smoothie: Efficient Inference of Spatial Co-expression Networks from Denoised Spatial Transcriptomics Data."
Chase Holdener and Iwijn De Vlaminck (2025). 
bioRxiv preprint: [https://doi.org/10.1101/2025.02.26.640406](https://doi.org/10.1101/2025.02.26.640406)

## Workflow
#### Smoothie has a three step algorithm:
1. Perform Gaussian smoothing on the spatial gene data to address noise and sparsity.
2. Calculate pairwise Pearson correlation coefficients (PCC) for each smoothed gene pair.
3. Generate a spatial co-expression network and identify spatial gene modules with graph-based clustering.

![Overview of Smoothie and Downstream Analyses](images/overview.png)

Smoothie can generate a gene pattern atlas for a provided spatial dataset.
We show the gene pattern atlas for the embryonic day 16.5 mouse embryo (Stereo-seq) below. (Data retreived from: https://db.cngb.org/stomics/mosta/)

![Smoothie's E16.5 Mouse Embryo Gene Pattern Atlas](images/E16_5_mouse_embryo_pattern_atlas.png)

## Requirements
#### Python(v3.10.12)
anndata(v0.9.2), igraph(v0.10.6), matplotlib(v3.7.1), numpy(v1.23.5), pandas(v1.5.3), scanpy(v1.9.3), scipy(v1.9.3).

The output network is visualized with Cytoscape.

## Installation
After installing required packages above, the Smoothie repository may be accessed with:
`git clone https://github.com/caholdener01/Smoothie.git`

## Usage
Example scripts for single dataset analysis and multi-dataset analysis are included in the folder `/examples`. 

### This repository is still under development.