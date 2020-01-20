library(monocle3)
packageVersion('monocle3')
library(reticulate)
import("louvain")

##### Read 10x Data ##############
# Provide the path to the Cell Ranger output.
cat <- load_cellranger_data("C:/Users/mengr/Dropbox/GGBC/cat/data_analysis", barcode_filtered = TRUE)
cat

cat = preprocess_cds(cat, method = "PCA", num_dim = 100)
plot_pc_variance_explained(cat)

# Umap-PCA
cat = reduce_dimension(cat, reduction_method="UMAP", preprocess_method = "PCA")

cat = cluster_cells(cat, resolution = 1e-3, verbose = T,
                    python_home = "C:/Users/theuser/location/.conda/envs/r-reticulate/python")
plot_cells(cat, reduction_method="UMAP")

cluster <- 