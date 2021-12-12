rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

u_mat1 <- bm[["rna.umap"]]@cell.embeddings
u_mat2 <- bm[["adt.umap"]]@cell.embeddings
tmp <- svd(t(u_mat1) %*% u_mat2)
rotation_mat <- tmp$u %*% t(tmp$v)
tmp <- u_mat2 %*% t(rotation_mat)
rownames(tmp) <- rownames(bm@meta.data)
colnames(tmp) <- colnames(bm[["adt.umap"]]@cell.embeddings)
bm[["adt.umap"]]@cell.embeddings <- tmp

n <- ncol(bm)
Seurat::DefaultAssay(bm) <- "RNA"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30)
bm <- Seurat::FindClusters(bm, resolution = 0.25)
tab_vec <- table(bm$RNA_snn_res.0.25)
round(tab_vec/n, 2)
rm_idx <- names(tab_vec[tab_vec/n < 0.02])
bm$RNA_snn_res.0.25[bm$RNA_snn_res.0.25 %in% rm_idx] <- NA

Seurat::DefaultAssay(bm) <- "ADT"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:18, reduction = "apca")
bm <- Seurat::FindClusters(bm, resolution = 0.25)
tab_vec <- table(bm$ADT_snn_res.0.25)
round(tab_vec/n, 2)
rm_idx <- names(tab_vec[tab_vec/n < 0.02])
bm$ADT_snn_res.0.25[bm$ADT_snn_res.0.25 %in% rm_idx] <- NA

plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                         group.by = "RNA_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq) RNA:\nMeta-clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_rna_metacluster.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "ADT_snn_res.0.25", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq) ADT:\nMeta-clustering"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_adt_metacluster.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

##############

color_vec2 <- rep("gray", nrow(bm@meta.data))
names(color_vec2) <- rownames(bm@meta.data)
plot1 <- Seurat::DimPlot(bm, 
                         reduction = "rna.umap", 
                         cols = color_vec2)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human Bone Marrow:\nCite-seq (RNA)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_rna_none.png"),
                plot1, device = "png", width = 4.5, height = 5, units = "in")


color_vec2 <- rep("gray", nrow(bm@meta.data))
names(color_vec2) <- rownames(bm@meta.data)
plot1 <- Seurat::DimPlot(bm, 
                         reduction = "adt.umap", 
                         cols = color_vec2)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human Bone Marrow:\nCite-seq (Protein)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1 <- plot1 + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_adt_none.png"),
                plot1, device = "png", width = 4.5, height = 5, units = "in")


####################################
load("../../../../out/Writeup14f/Writeup14f_citeseq_bm_dcca.RData")
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

vec1 <- dcca_res$score_1[,1]
vec2 <- dcca_res$score_2[,1]

color_df <- data.frame(celltype = sort(unique(bm@meta.data$celltype.l2)),
                       color = scales::hue_pal()(length(unique(bm@meta.data$celltype.l2))))
col_vec <- rep(NA, length(vec2))
uniq_celltypes <- sort(unique(bm@meta.data$celltype.l2))
for(celltype in uniq_celltypes){
  col_vec[which(bm@meta.data$celltype.l2 == celltype)] <- color_df[which(color_df$celltype == celltype),"color"]
}

png("../../../../out/figures/Writeup14g/example_scores.png",
    height = 1200, width = 600, units = "px", res = 300)
par(mar = c(0.5, 4, 0.5, 0.5))
plot(NA, ylim = range(c(vec1, vec2)), xlim = c(-.5, 1.5), 
     ylab = "Leading canonical score", xlab = "", bty = "n", xaxt = "n")
points(y = vec1, x = rep(0, length(vec1)) + runif(length(vec1), min = -.3, max = .3),
       pch = 16, col = col_vec, cex = 0.25)
points(y = vec2, x = rep(1, length(vec2)) + runif(length(vec1), min = -.3, max = .3), 
       pch = 16, col = col_vec, cex = 0.25)
graphics.off()
