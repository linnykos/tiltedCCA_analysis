rm(list=ls())
load("../../../../out/Writeup14i/Writeup14i_citeseq_bm_dcca.RData")

plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\nRNA cell-types"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_bm_rna_celltype.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")
plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\nADT cell-types"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_bm_adt_celltype.png"),
                plot1, device = "png", width = 6.5, height = 5, units = "in")


plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                         group.by = "RNA_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\nRNA large clusters"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_bm_rna_cluster.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")
plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "ADT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Cite-seq):\nADT large clusters"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14i/Writeup14i_citeseq_bm_adt_cluster.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#####################

color_df <- data.frame(celltype = sort(unique(bm$celltype.l2)),
                       color = scales::hue_pal()(length(unique(bm$celltype.l2))))
metacell_df <- as.data.frame(t(sapply(metacell_clustering, function(vec){
  celltype_tab <- table(bm$celltype.l2[vec])
  celltype_name <- names(celltype_tab)[which.max(celltype_tab)]
  col <- color_df[which(color_df$celltype == celltype_name), "color"]
  
  c(celltype = celltype_name, color = col)
})))

set.seed(10)
min_umap <- Seurat::RunUMAP(min_embedding)@cell.embeddings
png("../../../../out/figures/Writeup14i/Writeup14i_citeseq_bm_min_umap.png",
    height = 2000, width = 2000, res = 300, units = "px")
plot(min_umap[,1], min_umap[,2], col = metacell_df$color,
     asp = T, pch = 16, 
     xlab = "minUMAP_1", ylab = "minUMAP_2",
     main = "Human BM (Cite-seq):\nMin-embedding metacell averages")
graphics.off()

##########################

tmp <- form_kernel_matrix(mat = mat_1, K = 30, 
                          metacell_clustering = metacell_clustering, 
                          clustering_hierarchy = clustering_hierarchy_1,
                          symmetrize_func = "average")
dist_mat_1 <- tmp$dist_mat
kernel_mat_1 <- tmp$kernel_mat

tmp <- form_kernel_matrix(mat = mat_2, K = 18, 
                          metacell_clustering = metacell_clustering, 
                          clustering_hierarchy = clustering_hierarchy_2,
                          symmetrize_func = "average")
dist_mat_2 <- tmp$dist_mat
kernel_mat_2 <- tmp$kernel_mat

idx1 <- which(metacell_df$celltype == "CD4 Naive")
idx2 <- which(metacell_df$celltype == "CD8 Naive")

quantile(dist_mat_1[idx1, idx1]^2)
quantile(dist_mat_2[idx1, idx1]^2)

quantile(dist_mat_1[idx2, idx2]^2)
quantile(dist_mat_2[idx2, idx2]^2)

quantile(dist_mat_1[idx1, idx2]^2)
quantile(dist_mat_2[idx1, idx2]^2)

quantile(dist_mat_1[idx1, -c(idx1,idx2)]^2)
quantile(dist_mat_2[idx1, -c(idx1,idx2)]^2)

