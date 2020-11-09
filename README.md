# Single Cell data analysis using Seurat and Monocle3

Pipeline for the analysis of 10x single-cell RNA-sequencing data using Monocle3 and seurat. The pipeline including pesudo-time analysis and cell trajectory constructions. 

##  Raw data analysi using “Cellranger” and “10x Loupe Cell Browser”
Analysis pipeline using “Cellranger count” by 10x genomics

## Load dataset into R package “Monocle3” and preprocessing with PCA
The output of Cellranger analysis pipeline has three files: “barcodes”, “features” and “matrix”. The “barcodes” file contains all the cell names and the “feature” file contain the genes and gene short names. The “matrix” file includes the counts with the location of gene and cell. We load this data into Monocle3.

##  Cluster cells using UMAP and Louvain methods.
![Image of Yaktocat](https://github.com/zhanzmr/Single_cell_data_analysis/blob/master/figures/umap.PNG)

##  Find Marker gene for each cluster.
we run a test to find out what genes makes them different from one another. Then we can rank markers genes according to “pseudo R2” value. 

## Construct trajectories
![Image of Yaktocat](https://github.com/zhanzmr/Single_cell_data_analysis/blob/master/figures/trajectory.PNG)

## Oder cells in chronological order
we select cells that are in the first stage (Spermatogonial Stem Cells (SSC’s)-Stem Progenitor) and then construct the pseudo time analysis. By using the significant genes, we can find significant cells that contain significant genes. We will make a cutoff point here to subset the cells. The next step is to find the nearest nodes(points) on the trajectories for the significant cells. We order those nodes and construct the pseudo time analysis.

##  Pseudo time analysis

![Image of Yaktocat](https://github.com/zhanzmr/Single_cell_data_analysis/blob/master/figures/time.PNG)







