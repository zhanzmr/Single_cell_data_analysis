library(monocle3)
library(reticulate)
library(stringr)
library(Seurat)
library(ggplot2)
import("louvain")


##### Read 10x Data ##############
# Provide the path to the Cell Ranger output.
cat <- load_cellranger_data("C:/Users/mengr/Dropbox/GGBC/cat/data_analysis", barcode_filtered = TRUE)
cat

cat = preprocess_cds(cat, method = "PCA", num_dim = 100)
#plot_pc_variance_explained(cat)

#Umap-PCA
cat = reduce_dimension(cat, reduction_method="UMAP", preprocess_method = "PCA")

cat = cluster_cells(cat, resolution = 1e-4, verbose = T,
                    python_home = "C:/Users/mengr/location/.conda/envs/r-reticulate/python")
plot_cells(cat, reduction_method="UMAP")

#learn trajectory
cat <- learn_graph(cat)
plot_cells(cat, reduction_method="UMAP")


# Load the dataset
cat.data <- Read10X(data.dir = "C:/Users/mengr/Dropbox/GGBC/cat/data_analysis/outs/filtered_feature_bc_matrix")

#seurat normalization
cat_sc <- NormalizeData(cat.data, normalization.method = "LogNormalize", scale.factor = 10000)
cat_data <- as.data.frame(cat_sc)
cat_data$gene <- rownames(cat_data); rownames(cat_data) <- 1:nrow(cat_data)

# Read marker genes
marker <- read.csv("C:/Users/mengr/Dropbox/GGBC/cat/data_analysis/marker.csv", stringsAsFactors = F)
marker <- marker[,1:5]
cmlu <- str_split_fixed(marker$Genes, "- ", n = Inf)
cmlu[cmlu==""] <- NA
#cmlu <- na.omit(cmlu)

gene_names = cmlu[1,]
gene_names <- na.omit(gene_names)

## Subset matrix
sub_gene <- data.frame()
for (jj in 1:length(gene_names)){
  subset_data <- cat_data[cat_data$gene == gene_names[jj],]
#  subset_data <- subset_data[ , -which(names(subset_data) %in% c("gene"))]
  sub_gene <- rbind(sub_gene, subset_data)
}
#write.csv(sub_gene,"sub_gene.csv")
sub_gene <- sub_gene[ , -which(names(sub_gene) %in% c("gene"))]
t1 = unname(unlist(sub_gene))

t2 = t1[t1!=0]
x = seq(1,length(t2))

plot(x,t2)
boxplot(t2)

## Plot the matrix 
#x = seq(1,10000)
#gra <- ggplot(sub_gene, aes()) + 
#        geom_point(aes(x=seq, y=t()))

#########
string_cell_name <- c()
cutoff = 2.8

for (ii in 1:length(gene_names)){
  subset_data <- cat_data[cat_data$gene == gene_names[ii],]
  subset_data <- subset_data[ , -which(names(subset_data) %in% c("gene"))]
  
  temp <- colnames(subset_data[,subset_data >= cutoff])
  temp2 <- paste(temp, "-1",sep = "")
  string_cell_name <- c(string_cell_name,temp2)
  print(ii)
}

string = unique(string_cell_name)

get_earliest_principal_node <- function(cds,cell_ID){
  node_name = as.data.frame(igraph::V(principal_graph(cds)[["UMAP"]])$name)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.data.frame(closest_vertex)
  closest_vertex$cell <- rownames(closest_vertex); rownames(closest_vertex) <- 1:nrow(closest_vertex)
  names(closest_vertex) <- c("node","cell")
  index <- closest_vertex[closest_vertex$cell %in% cell_ID,]$node
  root_nodes <- as.character(node_name[index,])
  root_nodes
}


start <- get_earliest_principal_node(cat,string)
start <- unique(start)
cat <- order_cells(cat, root_pr_nodes=start)

plot_cells(cat,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5,
           cell_size = 0.8)

#cellid = as.numeric(rownames(dd[dd$cell %in% p,]))
#mm = closest_vertex[cell_ids,]
# cell_ids <- as.numeric(rownames(dd[dd$cell %in% temp2,]))
# closest_vertex <- cat@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
# closest_vertex <- as.matrix(closest_vertex[colnames(cat), ])
# root_pr_nodes <- igraph::V(principal_graph(cat)[["UMAP"]])$name[as.numeric(names(which.max(closest_vertex[cell_ids,])))]

# get_earliest_principal_node <- function(cds,cell_ID){
#   cell_ids <- as.numeric(rownames(dd[dd$cell %in% temp2,]))
#   closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
#                                                                              (which.max(table(closest_vertex[cell_ids,]))))]
#   root_pr_nodes
# }

# node_name = as.data.frame(igraph::V(principal_graph(cat)[["UMAP"]])$name)
# closest_vertex <- cat@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
# closest_vertex <- as.data.frame(closest_vertex)
# closest_vertex$cell <- rownames(closest_vertex); rownames(closest_vertex) <- 1:nrow(closest_vertex)
# names(closest_vertex) <- c("node","cell")
# index <- closest_vertex[closest_vertex$cell %in% temp2,]$node
# root_nodes <- as.character(node_name[index,])

#dd = as.data.frame(clusters(cat))
#dd$cell <- rownames(dd); rownames(dd) <- 1:nrow(dd)
#names(dd) = c("cluster", "cell")

#ID4 <- cat_data[cat_data$gene == "ID4",]
#ID4 <- ID4[ , -which(names(ID4) %in% c("gene"))]

#temp <- colnames(ID4[,ID4 >= 10])
#temp2 <- paste(temp, "-1",sep = "")

# row <- rownames(cat_data)
# col <- colnames(cat_data)
# # normalization-quantile
# cat_data <- as.matrix(cat.data)
# cat_data <- normalize.quantiles(cat_data, copy = TRUE)
# cat_data <- as.data.frame(cat_data)
# rownames(cat_data) <- row
# colnames(cat_data) <- col
#Deseq2

####################################################################################################################
####################################################################################################################

# ##### Read 10x Data ##############
# # Provide the path to the Cell Ranger output.
# cat <- load_cellranger_data("C:/Users/mengr/Dropbox/GGBC/cat/data_analysis", barcode_filtered = TRUE)
# cat
# 
# cat = preprocess_cds(cat, method = "PCA", num_dim = 100)
# #plot_pc_variance_explained(cat)
# 
# #Umap-PCA
# cat = reduce_dimension(cat, reduction_method="UMAP", preprocess_method = "PCA")
# 
# cat = cluster_cells(cat, resolution = 1e-4, verbose = T,
#                     python_home = "C:/Users/mengr/location/.conda/envs/r-reticulate/python")
# plot_cells(cat, reduction_method="UMAP")
# 
# #learn trajectory
# cat <- learn_graph(cat)
# plot_cells(cat, reduction_method="UMAP")
# 
# # Different Expression Test
# #Graph-autocorrelation analysis
# pr_graph_test_res = graph_test(cat, neighbor_graph="knn", cores=1)
# pr_deg_ids = row.names(subset(pr_graph_test_res, q_value < 0.05))
# 
# # Read marker genes
# marker <- read.csv("C:/Users/mengr/Dropbox/GGBC/cat/data_analysis/marker.csv", stringsAsFactors = F)
# marker <- marker[,1:5]
# cmlu <- str_split_fixed(marker$Genes, "- ", n = Inf)
# vec <- as.vector(t(cmlu))
# vec[vec==""] <- NA
# vec <- na.omit(vec)
# marker_t <- data.frame(order = c(rep(1,13),rep(2,1),rep(3,7),rep(4,17),rep(5,23),rep(6,26),rep(7,26),rep(8,27),rep(9,19),rep(10,7),
#                                  rep(11,12),rep(12,2)),gene = vec)
# 
# 
# mulc <- pr_graph_test_res[,5:7]
# rownames(mulc) <- 1:nrow(mulc)
# 
# clust <- as.data.frame(clusters(cat))
# clust$cell <- rownames(clust); rownames(clust) <- 1:nrow(clust); names(clust) = c("cluster","cell")

#plot gene





#cellid = as.numeric(rownames(dd[dd$cell %in% p,]))
#mm = closest_vertex[cell_ids,]
# cell_ids <- as.numeric(rownames(dd[dd$cell %in% temp2,]))
# closest_vertex <- cat@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
# closest_vertex <- as.matrix(closest_vertex[colnames(cat), ])
# root_pr_nodes <- igraph::V(principal_graph(cat)[["UMAP"]])$name[as.numeric(names(which.max(closest_vertex[cell_ids,])))]

# get_earliest_principal_node <- function(cds,cell_ID){
#   cell_ids <- as.numeric(rownames(dd[dd$cell %in% temp2,]))
#   closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
#   closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
#   root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
#                                                                              (which.max(table(closest_vertex[cell_ids,]))))]
#   root_pr_nodes
# }

# node_name = as.data.frame(igraph::V(principal_graph(cat)[["UMAP"]])$name)
# closest_vertex <- cat@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
# closest_vertex <- as.data.frame(closest_vertex)
# closest_vertex$cell <- rownames(closest_vertex); rownames(closest_vertex) <- 1:nrow(closest_vertex)
# names(closest_vertex) <- c("node","cell")
# index <- closest_vertex[closest_vertex$cell %in% temp2,]$node
# root_nodes <- as.character(node_name[index,])

#dd = as.data.frame(clusters(cat))
#dd$cell <- rownames(dd); rownames(dd) <- 1:nrow(dd)
#names(dd) = c("cluster", "cell")

#ID4 <- cat_data[cat_data$gene == "ID4",]
#ID4 <- ID4[ , -which(names(ID4) %in% c("gene"))]

#temp <- colnames(ID4[,ID4 >= 10])
#temp2 <- paste(temp, "-1",sep = "")

