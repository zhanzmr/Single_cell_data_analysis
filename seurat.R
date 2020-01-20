library(Seurat)

# Load the dataset
cat.data <- Read10X(data.dir = "C:/Users/mengr/Dropbox/GGBC/cat/data_analysis/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
cat_sc <- CreateSeuratObject(counts = cat.data, project = "cat_scRNA_seq", min.cells = 3, min.features = 200)
cat_sc

#cat_sc[["percent.mt"]] <- PercentageFeatureSet(cat_sc, pattern = "^MT-")
VlnPlot(cat_sc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

#cat_sc <- subset(cat_sc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
cat_sc <- NormalizeData(cat_sc, normalization.method = "LogNormalize", scale.factor = 10000)


cat_sc <- FindVariableFeatures(cat_sc, selection.method = "vst", nfeatures = 18000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cat_sc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cat_sc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(cat_sc)
cat_sc <- ScaleData(cat_sc, features = all.genes)


cat_sc <- RunPCA(cat_sc, features = VariableFeatures(object = cat_sc))
DimPlot(cat_sc, reduction = "pca")
DimHeatmap(cat_sc, dims = 1, cells = 500, balanced = TRUE)


cat_sc <- FindNeighbors(cat_sc, dims = 1:10)
cat_sc <- FindClusters(cat_sc, resolution = 0.2)

cat_sc <- RunUMAP(cat_sc, dims = 1:10)
# individual clusters
DimPlot(cat_sc, reduction = "umap", pt.size  = 1) + DarkTheme() + NoLegend()




