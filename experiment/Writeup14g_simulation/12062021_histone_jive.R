rm(list=ls())
load("../../out/simulation/Writeup14g_data.RData")

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

u_mat1 <- seurat_obj2[["umap1"]]@cell.embeddings
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
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/jive.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")
