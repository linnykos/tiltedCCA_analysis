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

###############################3

# let's start tinkering
pull_factor <- 0.9
push_factor <- 1.1
mat_2_new <- mat_2
seurat_obj2 <- seurat_obj

celltype_vec <- paste0("Blue ", 1:7)
mean_vec2 <- Matrix::colMeans(mat_2_new[which(!seurat_obj2$celltype_custom %in% celltype_vec),,drop = F])
for(celltype in celltype_vec){
  idx <- which(seurat_obj2$celltype_custom == celltype)
  mean_vec1 <- Matrix::colMeans(mat_2_new[idx,,drop = F])
  for(i in idx){
    mat_2_new[i,] <- mean_vec2 + push_factor*(mat_2_new[i,] - mean_vec2)
  }
  mean_vec1 <- Matrix::colMeans(mat_2_new[idx,,drop = F])
  for(i in idx){
    mat_2_new[i,] <- pull_factor*mean_vec1 + (1-pull_factor)*mat_2_new[i,]
  }
}

mat_2_new <- scale(mat_2_new, center = T, scale = F)
svd_res <- svd(mat_2_new)
svd_res$u <- svd_res$u[,1:20]
tmp1 <- median(svd_res$u[which(seurat_obj2$celltype_custom %in% celltype_vec),1])
tmp2 <- median(svd_res$u[which(!seurat_obj2$celltype_custom %in% celltype_vec),1])
sign_val <- ifelse(tmp1 < tmp2, -1, 1)
for(j in 1:length(celltype_vec)){
  idx <- which(seurat_obj2$celltype_custom == celltype_vec[j])
  svd_res$u[idx,1] <- svd_res$u[idx,1] + sign_val*0.002*j^(3/4)
}
svd_res$v <- svd_res$v[,1:20]
svd_res$d <- svd_res$d[1:20]
mat_2_new <- tcrossprod(svd_res$u %*% diag(svd_res$d), svd_res$v)
mat_2_new <- scale(mat_2_new, center = T, scale = F)
svd_res <- svd(mat_2_new)
svd_res$u <- svd_res$u[,1:20]
svd_res$v <- svd_res$v[,1:20]
svd_res$d <- svd_res$d[1:20]
svd_res$d[1] <- 27
svd_res$d[2:10] <- exp(seq(log(20), log(10), length.out = 9))
svd_res$d[11:20] <- exp(seq(log(10), log(1), length.out = 10))
mat_2_new <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
colnames(mat_2_new) <- colnames(mat_2)[1:ncol(mat_2_new)]
rownames(mat_2_new) <- rownames(mat_2)
seurat_obj2[["lsi"]] <- Seurat::CreateDimReducObject(mat_2_new, assay = "DNA")

##########################

mat_1_new <- mat_1
rowsum_vec <- apply(mat_1_new, 1, multiomicCCA:::.l2norm)
mat_1_new <- multiomicCCA:::.mult_vec_mat(1/rowsum_vec, mat_1_new)
pull_factor <- 0.7

idx <- which(seurat_obj2$celltype_l2 == "Cortical_neurons")
mean_vec <- Matrix::colMeans(mat_1_new[idx,,drop = F])
sd_val <- median(matrixStats::colSds(mat_1_new[idx,,drop = F]))/4
for(i in idx){
  mat_1_new[i,] <- pull_factor*mean_vec + (1-pull_factor)*mat_1_new[i,]
  mat_1_new[i,] <- mat_1_new[i,] + rnorm(ncol(mat_1_new), mean = 0, sd = sd_val)
}

mat_1_new <- scale(mat_1_new, center = T, scale = F)
svd_res <- svd(mat_1_new)
svd_res$u <- svd_res$u[,1:20]
svd_res$v <- svd_res$v[,1:20]
svd_res$d <- svd_res$d[1:20]
svd_res$d[1:10] <- seq(27, 15, length.out = 10)[1:10]
mat_1_new <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
colnames(mat_1_new) <- colnames(mat_2)[1:ncol(mat_1_new)]
rownames(mat_1_new) <- rownames(mat_2)
seurat_obj2[["pca"]] <- Seurat::CreateDimReducObject(mat_1_new, assay = "RNA")

###############################

# check the visualizations
set.seed(10)
seurat_obj2 <- Seurat::RunUMAP(seurat_obj2, dims = 1:20, reduction = "pca",
                               metric = "euclidean",
                               reduction.name = "umap1", 
                               reduction.key = "umapModality1_")

set.seed(10)
seurat_obj2 <- Seurat::RunUMAP(seurat_obj2, dims = 1:20, reduction = "lsi",
                               metric = "euclidean",
                               reduction.name = "umap2", 
                               reduction.key = "umapModality2_")

u_mat1 <- seurat_obj2[["umap1"]]@cell.embeddings
u_mat2 <- seurat_obj2[["umap2"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("umapModality2_1", "umapModality2_2")
seurat_obj2[["umap2"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "umap1", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Modality 1"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/umap1.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "umap2", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Modality 2"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/umap2.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

save(seurat_obj2, color_vec, color_df, 
     file = "../../out/simulation/Writeup14g_data.RData")
