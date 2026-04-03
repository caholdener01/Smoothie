"""
Smoothie

"A package for analyzing spatial transcriptomics data with Gaussian smoothing, spatial gene 
correlation network clustering, gene module visualization, and multi-sample integration."
"""

__version__ = "1.0.2"
__author__ = "Chase Holdener"
__license__ = "MIT"

# Import main functions from modules
from .gaussian_smoothing import run_parallelized_smoothing
from .network_analysis import make_spatial_network, make_geneset_spatial_network
from .spatial_correlation import compute_correlation_matrix, get_correlations_to_GOI
from .multi_sample_integration import (
    concatenate_smoothed_matrices, 
    run_second_order_correlation_analysis,
    plot_dataset_pair_stabilities,
    plot_gene_stability_distribution,
    plot_gene_stability,
    plot_module_stability
)
from .choosing_hyperparameters import select_clustering_params
from .plotting import (
    rotate_spatial,
    plot_gene,
    plot_modules,
    plot_gene_multisample,
    plot_modules_multisample
)
from .create_anndata_from_transcripts import create_anndata_from_transcripts
from .utils import suppress_warnings, enable_warnings, quiet_mode
from .shuffle_analysis import compute_shuffled_correlation_percentiles

__all__ = [
    # Smoothing functions
    "run_parallelized_smoothing",
    
    # Network analysis functions
    "make_spatial_network",
    "make_geneset_spatial_network",
    
    # Correlation functions
    "compute_correlation_matrix",
    "get_correlations_to_GOI",
    
    # Multi-sample integration functions
    "concatenate_smoothed_matrices",
    "run_second_order_correlation_analysis",
    "plot_dataset_pair_stabilities",
    "plot_gene_stability_distribution",
    "plot_gene_stability",
    "plot_module_stability",
    
    # Hyperparameter selection functions
    "select_clustering_params",
    
    # Plotting functions
    "rotate_spatial",
    "plot_gene",
    "plot_modules",
    "plot_gene_multisample",
    "plot_modules_multisample",
    
    # Imaging data conversion
    "create_anndata_from_transcripts",
    
    # Utilities
    "suppress_warnings",
    "enable_warnings",
    "quiet_mode",
    
    # Shuffle analysis (null distribution)
    "compute_shuffled_correlation_percentiles"
]
