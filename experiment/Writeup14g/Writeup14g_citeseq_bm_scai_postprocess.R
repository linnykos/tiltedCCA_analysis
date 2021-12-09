rm(list=ls())
load("../../../../out/Writeup14g/Writeup14g_citeseq_bm25_scai.RData")

for(j in 1:ncol(scai_res$H)){
  val <- max(scai_res$H[,j])
  scai_res$H[,j] <- scai_res$H[,j]/val
  scai_res$W1[,j] <- scai_res$W1[,j]*val
  scai_res$W2[,j] <- scai_res$W2[,j]*val
}

set.seed(10)
scai_umap <- Seurat::RunUMAP(scai_res$H, 
                             metric = "cosine",
                             reduction.key = "umapscAI_")
rownames(scai_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["scai_umap"]] <- Seurat::CreateDimReducObject(scai_umap@cell.embeddings)

u_mat1 <- bm[["rna.umap"]]@cell.embeddings
u_mat2 <- bm[["scai_umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(bm@meta.data)
colnames(tmp) <- colnames(bm[["scai_umap"]]@cell.embeddings)
bm[["scai_umap"]]@cell.embeddings <- tmp

plot1 <- Seurat::DimPlot(bm, reduction = "scai_umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\nscAI"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_scai_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")
