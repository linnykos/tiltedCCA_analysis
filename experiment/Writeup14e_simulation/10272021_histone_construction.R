rm(list=ls())
library(Seurat)
library(Signac)

# meant to be run on desktop
load("../../out/simulation/Writeup14e_simulation/10272021_H3K27me3_original.RData")

cell_list <- list(Non_neurons = c("Astro_Myoc", "Astro_Nnat", "OPC", "Microglia", 
                                  "Oligo_MFOL", "Oligo_MOL", "Endothelial", "Ependymal"),
                  Inhibitory_neurons = c("Pvalb", "Sst", "CGE"),
                  Cortical_neurons = c("CT", "NP", "L6", "L5", "L4", "L23", "PT"),
                  Hippocampal_neurons = c("DG", "Subiculum", "CA1", "CA23"))

mat_1b <- mat_1; colnames(mat_1b) <- paste0("gene-", 1:ncol(mat_1b))
seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1b))
mat_2b <- mat_2; colnames(mat_2b) <- paste0("dna-", 1:ncol(mat_2b))
seurat_obj[["DNA"]] <- Seurat::CreateAssayObject(counts = t(mat_2b))
seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(mat_1, assay = "RNA")
seurat_obj[["lsi"]] <- Seurat::CreateDimReducObject(mat_2, assay = "DNA")
seurat_obj[["celltype"]] <- celltype[,1]
celltype_l2 <- celltype[,1]
for(i in 1:length(cell_list)){
  celltype_l2[celltype_l2 %in% cell_list[[i]]] <- names(cell_list)[i]
}
seurat_obj[["celltype_l2"]] <- celltype_l2

uniq_celltypes <- sort(unique(seurat_obj@meta.data$celltype))
color_vec <- sapply(uniq_celltypes, function(i){
  color_df[which(color_df$celltype == i),"color"]
})


set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:50, reduction = "pca",
                              reduction.name = "umap")
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:25, reduction = "lsi",
                              metric = "euclidean",
                              reduction.name = "umap.dna")



Seurat::DimPlot(seurat_obj, reduction = "umap", cols = color_vec, group.by = "celltype")
Seurat::DimPlot(seurat_obj, reduction = "umap.dna", cols = color_vec, group.by = "celltype")

plot_scores_heatmap.list(list(mat_1[,1:25]), membership_vec = as.factor(seurat_obj[["celltype"]][,1]),
                         log_scale = T, scaling_power = 2)
plot_scores_heatmap.list(list(mat_2[,1:25]), membership_vec = as.factor(seurat_obj[["celltype_l2"]][,1]),
                         log_scale = T, scaling_power = 2)

##########################
##########################
##########################

# let's start tinkering
pull_factor <- 0.9
mat_2_new <- mat_2
seurat_obj2 <- seurat_obj

for(celltype in c("CT", "NP", "L6", "L5", "L4", "L23", "PT")){
  idx <- which(seurat_obj2$celltype == celltype)
  mean_vec <- Matrix::colMeans(mat_2_new[idx,,drop = F])
  for(i in idx){
    mat_2_new[i,] <- pull_factor*mean_vec + (1-pull_factor)*mat_2_new[i,]
  }
}

mat_2_new <- scale(mat_2_new, center = T, scale = F)
svd_res <- svd(mat_2_new)
svd_res$d[1:10] <- exp(seq(log(27), log(10), length.out = 10))
svd_res$d[11:50] <- exp(seq(log(10), log(1), length.out = 40))
mat_2_new <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
colnames(mat_2_new) <- colnames(mat_2)[1:ncol(mat_2_new)]
rownames(mat_2_new) <- rownames(mat_2)

seurat_obj2[["lsi"]] <- Seurat::CreateDimReducObject(mat_2_new, assay = "DNA")
set.seed(10)
seurat_obj2 <- Seurat::RunUMAP(seurat_obj2, dims = 1:20, reduction = "lsi",
                              metric = "euclidean",
                              reduction.name = "umap.dna",
                              reduction.key = "umapDNA_")
Seurat::DimPlot(seurat_obj2, reduction = "umap.dna", 
                cols = color_vec, 
                group.by = "celltype")

# plot_scores_heatmap.list(list(mat_2_new[,1:20]), membership_vec = as.factor(seurat_obj[["celltype_l2"]][,1]),
#                          log_scale = T, scaling_power = 2)

#######

mat_1_new <- mat_1
rowsum_vec <- apply(mat_1_new, 1, multiomicCCA:::.l2norm)
mat_1_new <- multiomicCCA:::.mult_vec_mat(1/rowsum_vec, mat_1_new)
pull_factor <- 0.7

set.seed(10)
idx <- which(seurat_obj2$celltype_l2 == "Cortical_neurons")
mean_vec <- Matrix::colMeans(mat_1_new[idx,,drop = F])
sd_val <- median(matrixStats::colSds(mat_1_new[idx,,drop = F]))/4
for(i in idx){
  mat_1_new[i,] <- pull_factor*mean_vec + (1-pull_factor)*mat_1_new[i,]
  mat_1_new[i,] <- mat_1_new[i,] + rnorm(ncol(mat_1_new), mean = 0, sd = sd_val)
}

mat_1_new <- scale(mat_1_new, center = T, scale = F)
svd_res <- svd(mat_1_new)
svd_res$d[1:10] <- seq(27, 13.5, by = -0.5)[1:10]
mat_1_new <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
colnames(mat_1_new) <- colnames(mat_2)[1:ncol(mat_1_new)]
rownames(mat_1_new) <- rownames(mat_2)

seurat_obj2[["pca"]] <- Seurat::CreateDimReducObject(mat_1_new, assay = "RNA")
set.seed(10)
seurat_obj2 <- Seurat::RunUMAP(seurat_obj2, dims = 1:20, reduction = "pca",
                               metric = "euclidean",
                               reduction.name = "umap", 
                               reduction.key = "umapRNA_")
Seurat::DimPlot(seurat_obj2, reduction = "umap", 
                cols = color_vec, 
                group.by = "celltype",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)

# plot_scores_heatmap.list(list(mat_1_new[,1:20]), membership_vec = as.factor(seurat_obj$celltype),
#                          log_scale = T, scaling_power = 2)

#####################

set.seed(10)
seurat_obj2 <- Seurat::FindMultiModalNeighbors(
  seurat_obj2, reduction.list = list("pca", "lsi"), 
  dims.list = list(1:20, 1:20), modality.weight.name = "RNA.weight"
)
set.seed(10)
seurat_obj2 <- Seurat::RunUMAP(seurat_obj2, 
                               nn.name = "weighted.nn", 
                               reduction.name = "wnn.umap", 
                               reduction.key = "wnnUMAP_")

Seurat::DimPlot(seurat_obj2, reduction = "wnn.umap", 
                cols = color_vec, 
                group.by = "celltype",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)

###########

n <- ncol(seurat_obj2)
svd_1 <- svd(seurat_obj2[["pca"]]@cell.embeddings)
embedding_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
svd_2 <- svd(seurat_obj2[["lsi"]]@cell.embeddings)
embedding_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
consensus_mat <- cbind(embedding_1, embedding_2)
rownames(consensus_mat) <- colnames(seurat_obj2)

set.seed(10)
consensus_umap <- Seurat::RunUMAP(consensus_mat, 
                                  metric = "euclidean",
                                  reduction.key = "consensus_")
seurat_obj2[["consensus.umap"]] <- Seurat::CreateDimReducObject(consensus_umap@cell.embeddings)

Seurat::DimPlot(seurat_obj2, reduction = "consensus.umap", 
                cols = color_vec, 
                group.by = "celltype",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)

# plot_scores_heatmap.list(list(consensus_mat), membership_vec = as.factor(seurat_obj$celltype),
#                          log_scale = T, scaling_power = 2)

##################################


source("../multiomicCCA_analysis/simulation/jive.R")
mat_1 <- seurat_obj2[["pca"]]@cell.embeddings
mat_2 <- seurat_obj2[["lsi"]]@cell.embeddings
jive_res <- jive(mat_1, mat_2, r = 20)

set.seed(10)
rownames(jive_res$embedding) <- colnames(seurat_obj2)
jive_umap <- Seurat::RunUMAP(jive_res$embedding, 
                                  metric = "euclidean",
                             reduction.key = "jive_")

seurat_obj2[["jive.umap"]] <- Seurat::CreateDimReducObject(jive_umap@cell.embeddings)

Seurat::DimPlot(seurat_obj2, reduction = "jive.umap", 
                cols = color_vec, 
                group.by = "celltype",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)

##########################
mat_1 <- seurat_obj2[["pca"]]@cell.embeddings
mat_2 <- seurat_obj2[["lsi"]]@cell.embeddings

snn_mat_1 <- .form_snn_mat(bool_intersect = T,
                         mat = mat_1, 
                         num_neigh = 30)
metacell_clustering_1 <- lapply(1:nrow(snn_mat_1), function(i){
  which(snn_mat_1[i,] != 0)
})
snn_mat_2 <- .form_snn_mat(bool_intersect = T,
                           mat = mat_2, 
                           num_neigh = 30)
metacell_clustering_2 <- lapply(1:nrow(snn_mat_2), function(i){
  which(snn_mat_2[i,] != 0)
})

dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      dims_1 = 1:20, dims_2 = 1:20,
                                      metacell_clustering_1 = metacell_clustering_1,
                                      metacell_clustering_2 = metacell_clustering_2,
                                      fix_tilt_perc = F)
dcca_res$cca_obj
dcca_res$df_percentage
dcca_res$tilt_perc
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, rank_c = 20)

set.seed(10)
dcca_common_umap <- Seurat::RunUMAP(cbind(dcca_decomp$common_mat_1, dcca_decomp$common_mat_2), 
                             metric = "euclidean",
                             reduction.key = "dccaCommon_")
seurat_obj2[["dcca_common"]] <- Seurat::CreateDimReducObject(dcca_common_umap@cell.embeddings)

Seurat::DimPlot(seurat_obj2, reduction = "dcca_common", 
                cols = color_vec, 
                group.by = "celltype",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)

set.seed(10)
dcca_distinct1_umap <- Seurat::RunUMAP(dcca_decomp$distinct_mat_1, 
                                    metric = "euclidean",
                                    reduction.key = "dccaDistinct1_")
seurat_obj2[["dcca_distinct1"]] <- Seurat::CreateDimReducObject(dcca_distinct1_umap@cell.embeddings)

Seurat::DimPlot(seurat_obj2, reduction = "dcca_distinct1", 
                cols = color_vec, 
                group.by = "celltype",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)

set.seed(10)
dcca_distinct2_umap <- Seurat::RunUMAP(dcca_decomp$distinct_mat_2, 
                                       metric = "euclidean",
                                       reduction.key = "dccaDistinct2_")
seurat_obj2[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(dcca_distinct2_umap@cell.embeddings)

Seurat::DimPlot(seurat_obj2, reduction = "dcca_distinct2", 
                cols = color_vec, 
                group.by = "celltype",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)

