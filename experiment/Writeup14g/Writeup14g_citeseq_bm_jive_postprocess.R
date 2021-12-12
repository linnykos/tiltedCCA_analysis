rm(list=ls())
load("../../../../out/Writeup14g/Writeup14g_citeseq_bm25_jive.RData")
dim(jive_res$embedding)

set.seed(10)
jive_umap <- Seurat::RunUMAP(jive_res$embedding, 
                             metric = "cosine",
                             reduction.key = "umapJive1_")
rownames(jive_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["jive_umap"]] <- Seurat::CreateDimReducObject(jive_umap@cell.embeddings)

u_mat1 <- bm[["rna.umap"]]@cell.embeddings
u_mat2 <- bm[["jive_umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(bm@meta.data)
colnames(tmp) <- colnames(bm[["jive_umap"]]@cell.embeddings)
bm[["jive_umap"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(bm, reduction = "jive_umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nJIVE"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_jive_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")
