rm(list=ls())
load("../../out/simulation/Writeup14g_data.RData")

source("../multiomicCCA_analysis/simulation/scai.R")
mat_1 <- seurat_obj2[["pca"]]@cell.embeddings 
for(j in 1:ncol(mat_1)){
  mat_1[,j] <- mat_1[,j] - min(mat_1[,j])
}
# svd_res <- svd(mat_1)
# for(j in 1:ncol(svd_res$u)){
#   svd_res$u[,j] <- svd_res$u[,j] + min(svd_res$u[,j])
#   svd_res$v[,j] <- svd_res$v[,j] + min(svd_res$v[,j])
# }
# tmp <- svd_res$u %*% diag(svd_res$d) %*% t(svd_res$v)
# rownames(tmp) <- rownames(mat_1)
# mat_1 <- tmp
mat_2 <- seurat_obj2[["lsi"]]@cell.embeddings
for(j in 1:ncol(mat_2)){
  mat_2[,j] <- mat_2[,j] - min(mat_2[,j])
}
# svd_res <- svd(mat_2)
# for(j in 1:ncol(svd_res$u)){
#   svd_res$u[,j] <- svd_res$u[,j] + min(svd_res$u[,j])
#   svd_res$v[,j] <- svd_res$v[,j] + min(svd_res$v[,j])
# }
# tmp <- svd_res$u %*% diag(svd_res$d) %*% t(svd_res$v)
# rownames(tmp) <- rownames(mat_2)
# mat_2 <- tmp
scai_res <- scai(mat_1, mat_2, r = 20, gamma = 0)
plot(scai_res$obj_vec)
quantile(scai_res$H)
quantile(scai_res$W1)
quantile(scai_res$W2)
par(mfrow = c(1,2))
image(t(mat_1))
image(t(scai_res2$H %*% t(scai_res2$W1)))
par(mfrow = c(1,2))
image(t(mat_2))
image(t(scai_res$H %*% t(scai_res$W2)))

apply(scai_res$H, 2, multiomicCCA:::.l2norm)
apply(scai_res$W1, 2, multiomicCCA:::.l2norm)
apply(scai_res$W2, 2, multiomicCCA:::.l2norm)
set.seed(10)
rownames(scai_res$H) <- colnames(seurat_obj2)
# tmp <- scai_res$H %*% t(scai_res$W2)
# svd_res <- svd(tmp); tmp <- svd_res$u %*% diag(svd_res$d)
scai_res2 <- scai_res
for(j in 1:ncol(scai_res2$H)){
  val <- max(scai_res2$H[,j])
  scai_res2$H[,j] <- scai_res2$H[,j]/val
  scai_res2$W1[,j] <- scai_res2$W1[,j]*val
  scai_res2$W2[,j] <- scai_res2$W2[,j]*val
}
tmp <- scai_res2$H
scai_umap <- Seurat::RunUMAP(tmp, 
                             metric = "euclidean",
                             reduction.key = "umapscai_")
seurat_obj2[["scai.umap"]] <- Seurat::CreateDimReducObject(scai_umap@cell.embeddings)

u_mat1 <- seurat_obj2[["umap1"]]@cell.embeddings
u_mat2 <- seurat_obj2[["scai.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(seurat_obj2@meta.data)
colnames(tmp) <- c("umapscai_1", "umapscai_2")
seurat_obj2[["scai.umap"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "scai.umap", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: scAI"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/scai.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

####################

umap1 <- Seurat::RunUMAP(mat_1, 
                         metric = "euclidean",
                         reduction.key = "umapscai_")
seurat_obj2[["umap1"]] <- Seurat::CreateDimReducObject(umap1@cell.embeddings)
plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "umap1", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Positive 1"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/umap1_positive.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


umap2 <- Seurat::RunUMAP(mat_2, 
                         metric = "euclidean",
                         reduction.key = "umapscai_")
seurat_obj2[["umap2"]] <- Seurat::CreateDimReducObject(umap2@cell.embeddings)
plot1 <- Seurat::DimPlot(seurat_obj2, reduction = "umap2", 
                         cols = color_vec, 
                         group.by = "celltype_custom",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: Positive 2"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/umap2_positive.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

