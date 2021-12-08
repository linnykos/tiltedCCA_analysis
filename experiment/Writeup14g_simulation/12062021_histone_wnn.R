rm(list=ls())
load("../../out/simulation/Writeup14g_data.RData")

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

u_mat1 <- seurat_obj2[["umap1"]]@cell.embeddings
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
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/wnn.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <-  Seurat::VlnPlot(seurat_obj2, features = "RNA.weight", 
                          cols = color_vec, 
                          group.by = "celltype_custom",
                          pt.size = 0.1) + Seurat::NoLegend()
plot1 <- plot1 + ggplot2::ggtitle(paste0("Pseudo-real data: WNN weights"))
ggplot2::ggsave(filename = paste0("../../out/figures/Writeup14g_simulation/wnn_weights.png"),
                plot1, device = "png", width = 8, height = 4, units = "in")
