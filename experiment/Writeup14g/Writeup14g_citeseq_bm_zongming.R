rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)
dims_1 <- 1:30; dims_2 <- 1:18; nn <- 30

max_mat <- zongming_embedding(mat_1, mat_2,
                              dims_1, dims_2, 
                              nn)

graph_obj <- SeuratObject::as.Graph(max_mat)

set.seed(10)
umap_res <- Seurat::RunUMAP(graph_obj, 
                            metric = "euclidean",
                            assay = "RNA",
                            reduction.key = "umapMax_")
rownames(umap_res@cell.embeddings) <- rownames(mat_1)
bm[["max.umap"]] <- Seurat::CreateDimReducObject(umap_res@cell.embeddings)

u_mat1 <- bm[["rna.umap"]]@cell.embeddings
u_mat2 <- bm[["max.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(bm@meta.data)
colnames(tmp) <- colnames(bm[["max.umap"]]@cell.embeddings)
bm[["max.umap"]]@cell.embeddings <- tmp


plot1 <- Seurat::DimPlot(bm, reduction = "max.umap", 
                         group.by = "celltype.l2",
                         label = TRUE, 
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nMaximum distance"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_maxdist.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")



