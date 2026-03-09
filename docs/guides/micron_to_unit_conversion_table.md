For correct Gaussian smoothing of your spatial transcriptomics dataset, you must know the conversion factor between microns to spatial units. This table provides conversion factors for several popular spatial platforms. 

#### **Sequencing-based Spatial Platforms**

| Spatial Platform | Micron to Unit Conversion (1 micron = X units) | Notes |
| :--- | :--- | :--- |
| **Slide-seq (Curio-Seeker)** | 1.5456 | 1 micron = 1.5456 Slide-seq units |
| **Stereo-seq (STOmics)** | 2 | 1 micron = 2 Stereo-seq units |
| **Visium-HD (10X Genomics)** | (1 / `microns_per_pixel`) | You can find the VisiumHD `microns_per_pixel` conversion value in scalefactors_json.json under the directory binned_outputs/square_XXX/spatial. This varies per experiment. |

#### **Imaging-based Spatial Platforms**

| Spatial Platform | Micron to Unit Conversion (1 micron = X units) | Notes |
| :--- | :--- | :--- |
| **Xenium (10X Genomics)** | 1 | transcripts.parquet stores transcript coordinates in microns. |
| **MERFISH (MERSCOPE / MERSCOPE Ultra)** | 1 | detected_transcripts.parquet stores transcript coordinates in microns. |
| **CosMx SMI (Bruker)** | 8.31393 *or* 0.001 | Units may be stored as pixels (_px) or millimeters (_mm), indicated in the column names. (1 micron = 8.31393 pixels) |
| **seqFISH (Spatial Genomics)** | Varies per experiment | Check image information to get this value. |

!!! warning
    Be sure to double check the conversion value is correct for your dataset! 
    
If you wish to add new conversion values for spatial omics platforms not listed in this table, please feel free to add a github issue.