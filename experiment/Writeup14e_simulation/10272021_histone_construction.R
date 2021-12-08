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
color_name_vec <- c("Orange", "Purple", "Blue", "Green")
celltype_custom <- celltype[,1]
celltype_mapping <- color_df$celltype
for(i in 1:length(cell_list)){
  for(j in 1:length(cell_list[[i]])){
    celltype_custom[which(celltype_custom == cell_list[[i]][j])] <- paste0(color_name_vec[i], " ", j)
    celltype_mapping[which(celltype_mapping == cell_list[[i]][j])] <- paste0(color_name_vec[i], " ", j)
  }
}
color_df$custom <- celltype_mapping

mat_1b <- mat_1; colnames(mat_1b) <- paste0("gene-", 1:ncol(mat_1b))
seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1b))
mat_2b <- mat_2; colnames(mat_2b) <- paste0("dna-", 1:ncol(mat_2b))
seurat_obj[["DNA"]] <- Seurat::CreateAssayObject(counts = t(mat_2b))
seurat_obj[["pca"]] <- Seurat::CreateDimReducObject(mat_1, assay = "RNA")
seurat_obj[["lsi"]] <- Seurat::CreateDimReducObject(mat_2, assay = "DNA")
seurat_obj[["celltype"]] <- celltype[,1]
seurat_obj[["celltype_custom"]] <- celltype_custom
celltype_l2 <- celltype[,1]
for(i in 1:length(cell_list)){
  celltype_l2[celltype_l2 %in% cell_list[[i]]] <- names(cell_list)[i]
}
seurat_obj[["celltype_l2"]] <- celltype_l2

uniq_celltypes <- sort(unique(seurat_obj@meta.data$celltype_custom))
color_vec <- sapply(uniq_celltypes, function(i){
  color_df[which(color_df$custom == i),"color"]
})


set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:50, reduction = "pca",
                              reduction.name = "umap")
set.seed(10)
seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:25, reduction = "lsi",
                              metric = "euclidean",
                              reduction.name = "umap.dna")

Seurat::DimPlot(seurat_obj, reduction = "umap", cols = color_vec, group.by = "celltype_custom")
Seurat::DimPlot(seurat_obj, reduction = "umap.dna", cols = color_vec, group.by = "celltype_custom")

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
svd_res$d[1] <- 27
svd_res$d[2:10] <- exp(seq(log(20), log(10), length.out = 9))
svd_res$d[11:50] <- exp(seq(log(10), log(1), length.out = 40))
mat_2_new <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
colnames(mat_2_new) <- colnames(mat_2)[1:ncol(mat_2_new)]
rownames(mat_2_new) <- rownames(mat_2)

seurat_obj2[["lsi"]] <- Seurat::CreateDimReducObject(mat_2_new, assay = "DNA")
set.seed(10)
seurat_obj2 <- Seurat::RunUMAP(seurat_obj2, dims = 1:20, reduction = "lsi",
                              metric = "euclidean",
                              reduction.name = "umap.dna",
                              reduction.key = "umapModality2_")
set.seed(10)
Seurat::DefaultAssay(seurat_obj2) <- "DNA"
seurat_obj2 <- Seurat::FindNeighbors(seurat_obj2, reduction = "lsi",
                                     assay = "DNA",
                                     dims = 1:20)
seurat_obj2 <- Seurat::FindClusters(seurat_obj2, resolution = 1)
table(seurat_obj2$DNA_snn_res.1)
vec <- seurat_obj2$DNA_snn_res.1
vec[which(vec == "3")] <- "2"
vec <- as.factor(as.numeric(vec))
levels(vec) <- 1:length(levels(vec))
seurat_obj2$DNA_snn_res.1 <- vec

# plot_scores_heatmap.list(list(mat_2_new[,1:20]), membership_vec = as.factor(seurat_obj[["celltype_l2"]][,1]),
#                          log_scale = T, scaling_power = 5)

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
svd_res$d[1:10] <- seq(27, 15, length.out = 10)[1:10]
mat_1_new <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
colnames(mat_1_new) <- colnames(mat_2)[1:ncol(mat_1_new)]
rownames(mat_1_new) <- rownames(mat_2)

seurat_obj2[["pca"]] <- Seurat::CreateDimReducObject(mat_1_new, assay = "RNA")
set.seed(10)
seurat_obj2 <- Seurat::RunUMAP(seurat_obj2, dims = 1:20, reduction = "pca",
                               metric = "euclidean",
                               reduction.name = "umap", 
                               reduction.key = "umapModality1_")

set.seed(10)
Seurat::DefaultAssay(seurat_obj2) <- "RNA"
seurat_obj2 <- Seurat::FindNeighbors(seurat_obj2, reduction = "pca",
                                     assay = "RNA",
                                     dims = 1:20)
seurat_obj2 <- Seurat::FindClusters(seurat_obj2, resolution = 1)
table(seurat_obj2$RNA_snn_res.1)
vec <- seurat_obj2$RNA_snn_res.1
vec <- as.character(vec)
vec[which(vec %in% c("1", "4"))] <- "0"
vec <- as.factor(as.numeric(vec))
levels(vec) <- 1:length(levels(vec))
seurat_obj2$RNA_snn_res.1 <- vec


# plot_scores_heatmap.list(list(mat_1_new[,1:20]), membership_vec = as.factor(seurat_obj$celltype),
#                          log_scale = T, scaling_power = 2)

############################

u_mat1 <- seurat_obj2[["umap"]]@cell.embeddings
u_mat2 <- seurat_obj2[["umap.dna"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("umapModality2_1", "umapModality2_2")
seurat_obj2[["umap.dna"]]@cell.embeddings <- tmp


##################################
##################################
##################################

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "umap.dna", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Modality 2"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/modality2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "umap.dna", 
                         group.by = "DNA_snn_res.1",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Modality 2"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/modality2_estimated.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")


color_vec2 <- rep("gray", length(color_vec))
names(color_vec2) <- names(color_vec)
plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "umap.dna", cols = color_vec2)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Modality 2"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/modality2_none.png"),
                plot1, device = "png", width = 4.5, height = 5, units = "in")


##################################

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "umap", 
                cols = color_vec, 
                group.by = "celltype_custom",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Modality 1"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/modality1.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "umap", 
                         group.by = "RNA_snn_res.1",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Modality 1"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/modality1_estimated.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

color_vec2 <- rep("gray", length(color_vec))
names(color_vec2) <- names(color_vec)
plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "umap", cols = color_vec2)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Modality 1"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/modality1_none.png"),
                plot1, device = "png", width = 4.5, height = 5, units = "in")


##################################
##################################
##################################

set.seed(10)
seurat_obj2 <- Seurat::FindMultiModalNeighbors(
  seurat_obj2, reduction.list = list("pca", "lsi"), 
  dims.list = list(1:20, 1:20), modality.weight.name = "RNA.weight"
)
set.seed(10)
seurat_obj2 <- Seurat::RunUMAP(seurat_obj2, 
                               nn.name = "weighted.nn", 
                               reduction.name = "wnn.umap", 
                               reduction.key = "umapWNN_")

u_mat1 <- seurat_obj2[["umap"]]@cell.embeddings
u_mat2 <- seurat_obj2[["wnn.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("umapWNN_1", "umapWNN_2")
seurat_obj2[["wnn.umap"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "wnn.umap", 
                cols = color_vec, 
                group.by = "celltype_custom",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: WNN"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/wnn.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <-  Seurat::VlnPlot(seurat_obj2, features = "RNA.weight", 
                          cols = color_vec, 
                          group.by = "celltype_custom",
                          pt.size = 0.1) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: WNN weights"))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/wnn_weights.png"),
                plot1, device = "png", width = 8, height = 5, units = "in")

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
                                  reduction.key = "umapConsensusPCA_")
seurat_obj2[["consensus.umap"]] <- Seurat::CreateDimReducObject(consensus_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap"]]@cell.embeddings
u_mat2 <- seurat_obj2[["consensus.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("umapConsensusPCA_1", "umapConsensusPCA_2")
seurat_obj2[["consensus.umap"]]@cell.embeddings <- tmp


plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "consensus.umap", 
                cols = color_vec, 
                group.by = "celltype_custom",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Consensus PCA"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/consensuspca.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


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
                             reduction.key = "umapJive_")

seurat_obj2[["jive.umap"]] <- Seurat::CreateDimReducObject(jive_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap"]]@cell.embeddings
u_mat2 <- seurat_obj2[["jive.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("umapJive_1", "umapJive_2")
seurat_obj2[["jive.umap"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "jive.umap", 
                cols = color_vec, 
                group.by = "celltype_custom",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: JIVE"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/jive.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


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
                                      fix_tilt_perc = 0.5)
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, rank_c = 20)

set.seed(10)
dcca_common_umap <- Seurat::RunUMAP(cbind(dcca_decomp$common_mat_1, dcca_decomp$common_mat_2), 
                                    metric = "euclidean",
                                    reduction.key = "dccaCommon_")
seurat_obj2[["dcca_common"]] <- Seurat::CreateDimReducObject(dcca_common_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap"]]@cell.embeddings
u_mat2 <- seurat_obj2[["dcca_common"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("dccaCommon_1", "dccaCommon_2")
seurat_obj2[["dcca_common"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "dcca_common", 
                cols = color_vec, 
                group.by = "celltype_custom",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: D-CCA (Common)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/dcca_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

color_vec2 <- rep("gray", length(color_vec))
names(color_vec2) <- names(color_vec)
plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "dcca_common", 
                         cols = color_vec2)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: D-CCA (Common)"))
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/dcca_common_none.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")


####################################

snn_mat_3 <- .form_snn_mat(bool_intersect = T,
                           mat = dcca_decomp$common_score, 
                           num_neigh = 30)
blue_idx <- which(seurat_obj2$celltype_custom == "Blue 1")
nn_size <- sapply(blue_idx, function(i){
  min(c(length(which(snn_mat_1[i,] != 0)), 
        length(which(snn_mat_2[i,] != 0)),
        length(which(snn_mat_3[i,] != 0))))
})
i <- blue_idx[which.max(nn_size)]

n <- nrow(mat_1)
vec <- rep("none", n)
vec[i] <- "target"
vec[which(snn_mat_1[i,] != 0)] <- "neigh"
vec[intersect(which(snn_mat_1[i,] != 0), which(snn_mat_3[i,] != 0))] <- "intersect"
seurat_obj2$mod1_neigh <- vec

vec <- rep("none", n)
vec[i] <- "target"
vec[which(snn_mat_2[i,] != 0)] <- "neigh"
vec[intersect(which(snn_mat_2[i,] != 0), which(snn_mat_3[i,] != 0))] <- "intersect"
seurat_obj2$mod2_neigh <- vec

vec <- rep("none", n)
vec[i] <- "target"
vec[which(snn_mat_3[i,] != 0)] <- "neighdcca"
vec[intersect(which(snn_mat_1[i,] != 0), which(snn_mat_3[i,] != 0))] <- "intersect"
seurat_obj2$dcca_neigh1 <- vec

vec <- rep("none", n)
vec[i] <- "target"
vec[which(snn_mat_3[i,] != 0)] <- "neighdcca"
vec[intersect(which(snn_mat_2[i,] != 0), which(snn_mat_3[i,] != 0))] <- "intersect"
seurat_obj2$dcca_neigh2 <- vec

length(intersect(which(snn_mat_1[i,] != 0), which(snn_mat_3[i,] != 0)))/length(unique(c(which(snn_mat_1[i,] != 0), which(snn_mat_3[i,] != 0))))
length(intersect(which(snn_mat_2[i,] != 0), which(snn_mat_3[i,] != 0)))/length(unique(c(which(snn_mat_2[i,] != 0), which(snn_mat_3[i,] != 0))))

color_vec3 <- c(rgb(127, 127, 127, max = 255), 
                rgb(212, 65, 86, max = 255), 
                rgb(72, 161, 217, max = 255), 
                rgb(249, 97, 255, max = 255), 
                rgb(77, 194, 61, max = 255))
names(color_vec3) <- c("none", "neigh", "neighdcca", "intersect", "target")

umap_mat <- seurat_obj2[["umap"]]@cell.embeddings
png("../../out/simulation/Writeup14e_simulation/modality1_none2.png", 
    height = 1650, width = 1500, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[,1], umap_mat[,2], col = color_vec3["none"], pch = 16, cex = 0.75)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/modality1_target.png", 
    height = 1650, width = 1500, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/modality1_neigh_original.png", 
    height = 1650, width = 1500, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$mod1_neigh == "neigh")
idx2 <- which(seurat_obj2$mod1_neigh == "intersect")
points(umap_mat[c(idx1,idx2),1], umap_mat[c(idx1,idx2),2], col = "white", pch = 16, cex = 3)
points(umap_mat[c(idx1,idx2),1], umap_mat[c(idx1,idx2),2], col = color_vec3["neigh"], pch = 16, cex = 2)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/modality1_neigh.png", 
    height = 1650, width = 1500, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$mod1_neigh == "neigh")
idx2 <- which(seurat_obj2$mod1_neigh == "intersect")
points(umap_mat[idx1,1], umap_mat[idx1,2], col = "gray", pch = 16, cex = 3)
points(umap_mat[idx1,1], umap_mat[idx1,2], col = color_vec3["neigh"], pch = 16, cex = 2)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = "white", pch = 16, cex = 3)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = color_vec3["intersect"], pch = 16, cex = 2)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/modality1_neighonly.png", 
    height = 1650, width = 1500, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$mod1_neigh == "neigh")
idx2 <- which(seurat_obj2$mod1_neigh == "intersect")
points(umap_mat[c(idx1,idx2,i),1], umap_mat[c(idx1,idx2,i),2], col = "white", pch = 16, cex = 3)
points(umap_mat[c(idx1,idx2,i),1], umap_mat[c(idx1,idx2,i),2], col = color_vec3["neigh"], pch = 16, cex = 2)
graphics.off()
tmp <- apply(umap_mat[c(idx1,idx2),], 2, range)
xlim <- tmp[,1]; ylim <- tmp[,2]
png("../../out/simulation/Writeup14e_simulation/modality1_neigh_zoom_original.png", 
    height = 600, width = 600, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = xlim, ylim = ylim,
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$mod1_neigh == "neigh")
idx2 <- which(seurat_obj2$mod1_neigh == "intersect")
points(umap_mat[c(idx1,idx2),1], umap_mat[c(idx1,idx2),2], col = "white", pch = 16, cex = 3)
points(umap_mat[c(idx1,idx2),1], umap_mat[c(idx1,idx2),2], col = color_vec3["neigh"], pch = 16, cex = 2)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/modality1_neigh_zoom.png", 
    height = 600, width = 600, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = xlim, ylim = ylim,
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$mod1_neigh == "neigh")
idx2 <- which(seurat_obj2$mod1_neigh == "intersect")
points(umap_mat[idx1,1], umap_mat[idx1,2], col = "gray", pch = 16, cex = 3)
points(umap_mat[idx1,1], umap_mat[idx1,2], col = color_vec3["neigh"], pch = 16, cex = 2)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = "white", pch = 16, cex = 3)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = color_vec3["intersect"], pch = 16, cex = 2)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()



umap_mat <- seurat_obj2[["umap.dna"]]@cell.embeddings
png("../../out/simulation/Writeup14e_simulation/modality2_none2.png", 
    height = 1650, width = 1500, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[,1], umap_mat[,2], col = color_vec3["none"], pch = 16, cex = 0.75)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/modality2_target.png", 
    height = 1650, width = 1500, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/modality2_neigh.png", 
    height = 1650, width = 1500, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$mod2_neigh == "neigh")
idx2 <- which(seurat_obj2$mod2_neigh == "intersect")
points(umap_mat[idx1,1], umap_mat[idx1,2], col = "gray", pch = 16, cex = 3)
points(umap_mat[idx1,1], umap_mat[idx1,2], col = color_vec3["neigh"], pch = 16, cex = 2)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = "white", pch = 16, cex = 3)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = color_vec3["intersect"], pch = 16, cex = 2)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/modality2_neighonly.png", 
    height = 1650, width = 1500, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$mod1_neigh == "neigh")
idx2 <- which(seurat_obj2$mod1_neigh == "intersect")
points(umap_mat[c(idx1,idx2,i),1], umap_mat[c(idx1,idx2,i),2], col = "white", pch = 16, cex = 3)
points(umap_mat[c(idx1,idx2,i),1], umap_mat[c(idx1,idx2,i),2], col = color_vec3["neigh"], pch = 16, cex = 2)
graphics.off()
tmp <- apply(umap_mat[c(idx1,idx2),], 2, range)
xlim <- tmp[,1]; ylim <- tmp[,2]
png("../../out/simulation/Writeup14e_simulation/modality2_neigh_zoom.png", 
    height = 600, width = 600, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = xlim, ylim = ylim,
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$mod2_neigh == "neigh")
idx2 <- which(seurat_obj2$mod2_neigh == "intersect")
points(umap_mat[idx1,1], umap_mat[idx1,2], col = "gray", pch = 16, cex = 3)
points(umap_mat[idx1,1], umap_mat[idx1,2], col = color_vec3["neigh"], pch = 16, cex = 2)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = "white", pch = 16, cex = 3)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = color_vec3["intersect"], pch = 16, cex = 2)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()


umap_mat <- seurat_obj2[["dcca_common"]]@cell.embeddings
png("../../out/simulation/Writeup14e_simulation/dcca_common_none.png", 
    height = 1500, width = 1600, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[,1], umap_mat[,2], col = color_vec3["none"], pch = 16, cex = 0.75)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/dcca_common_target.png", 
    height = 1500, width = 1600, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/dcca_common_neighonly.png", 
    height = 1500, width = 1600, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$dcca_neigh1 == "neighdcca")
idx2 <- which(seurat_obj2$dcca_neigh1 == "intersect")
points(umap_mat[c(idx1,idx2),1], umap_mat[c(idx1,idx2),2], col = "white", pch = 16, cex = 3)
points(umap_mat[c(idx1,idx2),1], umap_mat[c(idx1,idx2),2], col = color_vec3["neighdcca"], pch = 16, cex = 2)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/dcca_common_neigh1.png", 
    height = 1500, width = 1600, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$dcca_neigh1 == "neighdcca")
idx2 <- which(seurat_obj2$dcca_neigh1 == "intersect")
points(umap_mat[idx1,1], umap_mat[idx1,2], col = "gray", pch = 16, cex = 3)
points(umap_mat[idx1,1], umap_mat[idx1,2], col = color_vec3["neighdcca"], pch = 16, cex = 2)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = "white", pch = 16, cex = 3)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = color_vec3["intersect"], pch = 16, cex = 2)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()
png("../../out/simulation/Writeup14e_simulation/dcca_common_neigh2.png", 
    height = 1500, width = 1600, units = "px", res = 300)
par(mar = c(0,0,0,0))
plot(NA, xlim = range(umap_mat[,1]),
     ylim = range(umap_mat[,2]),
     xaxt = "n", yaxt = "n", bty = "n")
points(umap_mat[-i,1], umap_mat[-i,2], col = color_vec3["none"], pch = 16, cex = 0.75)
idx1 <- which(seurat_obj2$dcca_neigh2 == "neighdcca")
idx2 <- which(seurat_obj2$dcca_neigh2 == "intersect")
points(umap_mat[idx1,1], umap_mat[idx1,2], col = "gray", pch = 16, cex = 3)
points(umap_mat[idx1,1], umap_mat[idx1,2], col = color_vec3["neighdcca"], pch = 16, cex = 2)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = "white", pch = 16, cex = 3)
points(umap_mat[idx2,1], umap_mat[idx2,2], col = color_vec3["intersect"], pch = 16, cex = 2)
points(umap_mat[i,1], umap_mat[i,2], col = "white", pch = 16, cex = 3)
points(umap_mat[i,1], umap_mat[i,2], col = color_vec3["target"], pch = 16, cex = 2)
graphics.off()

####################################
####################################
####################################

dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      dims_1 = 1:20, dims_2 = 1:20,
                                      metacell_clustering_1 = metacell_clustering_1,
                                      metacell_clustering_2 = metacell_clustering_2,
                                      fix_tilt_perc = F)
dcca_res$tilt_perc
dcca_res$df_percentage

dcca_res2 <- fine_tuning(dcca_res, verbose = T)
dcca_res2$cca_obj
dcca_res2$tilt_perc
svd(dcca_res$common_score)$d[1]
svd(dcca_res2$common_score - dcca_res$common_score)$d[1]
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res2, rank_c = 20)

set.seed(10)
dcca_common_umap <- Seurat::RunUMAP(cbind(dcca_decomp$common_mat_1, dcca_decomp$common_mat_2), 
                             metric = "euclidean",
                             reduction.key = "dccaCommon_")
seurat_obj2[["dcca_common"]] <- Seurat::CreateDimReducObject(dcca_common_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap"]]@cell.embeddings
u_mat2 <- seurat_obj2[["dcca_common"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("dccaCommon_1", "dccaCommon_2")
seurat_obj2[["dcca_common"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "dcca_common", 
                cols = color_vec, 
                group.by = "celltype_custom",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Tilted-CCA\n(Common)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/tiltedcca_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


set.seed(10)
dcca_distinct1_umap <- Seurat::RunUMAP(dcca_decomp$distinct_mat_1, 
                                    metric = "euclidean",
                                    reduction.key = "dccaDistinct1_")
seurat_obj2[["dcca_distinct1"]] <- Seurat::CreateDimReducObject(dcca_distinct1_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap"]]@cell.embeddings
u_mat2 <- seurat_obj2[["dcca_distinct1"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("dccaDistinct1_1", "dccaDistinct1_2")
seurat_obj2[["dcca_distinct1"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "dcca_distinct1", 
                cols = color_vec, 
                group.by = "celltype_custom",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Tilted-CCA\n(Distinct 1)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/tiltedcca_distinct1.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


set.seed(10)
dcca_distinct2_umap <- Seurat::RunUMAP(dcca_decomp$distinct_mat_2, 
                                       metric = "euclidean",
                                       reduction.key = "dccaDistinct2_")
seurat_obj2[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(dcca_distinct2_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap"]]@cell.embeddings
u_mat2 <- seurat_obj2[["dcca_distinct2"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("dccaDistinct2_1", "dccaDistinct2_2")
seurat_obj2[["dcca_distinct2"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "dcca_distinct2", 
                cols = color_vec, 
                group.by = "celltype_custom",
                label = TRUE, 
                repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Tilted-CCA\n(Distinct 2)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/simulation/Writeup14e_simulation/tiltedcca_distinct2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

##################################


vec1 <- dcca_res$score_1[,1]
vec2 <- dcca_res$score_2[,1]
col_vec <- rep(NA, length(vec1))
uniq_celltypes <- sort(unique(seurat_obj2@meta.data$celltype))
for(celltype in uniq_celltypes){
  col_vec[which(seurat_obj2@meta.data$celltype == celltype)] <- color_df[which(color_df$celltype == celltype),"color"]
}

png("../../out/simulation/Writeup14e_simulation/example_scores.png",
    height = 1200, width = 600, units = "px", res = 300)
par(mar = c(0.5, 4, 0.5, 0.5))
plot(NA, ylim = range(c(vec1, vec2)), xlim = c(-.5, 1.5), 
     ylab = "Leading canonical score", xlab = "", bty = "n", xaxt = "n")
points(y = vec1, x = rep(0, length(vec1)) + runif(length(vec1), min = -.3, max = .3),
       pch = 16, col = col_vec, cex = 0.25)
points(y = vec2, x = rep(1, length(vec2)) + runif(length(vec1), min = -.3, max = .3), 
       pch = 16, col = col_vec, cex = 0.25)
graphics.off()


###################################

n <- nrow(mat_1)
svd_1 <- svd(mat_1); head(svd_1$d)
svd_2 <- svd(mat_2); head(svd_2$d)
resid_1 <- (diag(n) - tcrossprod(svd_2$u)) %*% mat_1
resid_svd_1 <- svd(resid_1); head(resid_svd_1$d)
resid_2 <- (diag(n) - tcrossprod(svd_1$u)) %*% mat_2
resid_svd_2 <- svd(resid_2); head(resid_svd_2$d)
