rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/Writeup14d/Writeup14d_citeseq_pbmc228_dcca.RData")

keep_vec <- rep(0, ncol(pbmc))
keep_vec[which(rownames(pbmc@meta.data) %in% rownames(rna_embeddings[[1]]))] <- 1
pbmc[["keep"]] <- keep_vec
pbmc <- subset(pbmc, keep == 1)

pbmc[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "common_", assay = "SCT")
pbmc[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "distinct_", assay = "SCT")
pbmc[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "everything_", assay = "SCT")
rm(list = c("c_eig", "d_eig", "e_eig", "rna_embeddings", "rna_frnn")); gc(T)

pbmc[["common2"]] <- Seurat::CreateDimReducObject(embedding = protein_embeddings[[1]], key = "common2_", assay = "SCT")
pbmc[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = protein_embeddings[[2]], key = "distinct2_", assay = "SCT")
pbmc[["everything2"]] <- Seurat::CreateDimReducObject(embedding = protein_embeddings[[3]], key = "everything2_", assay = "SCT")
rm(list = c("c_eig2", "d_eig2", "e_eig2", "protein_embeddings", "protein_frnn")); gc(T)

pbmc[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, key = "combined_", assay = "SCT")
rm(list = c("combined_common_umap")); gc(T)

embedding_1 <- pbmc[["pca"]]@cell.embeddings
embedding_2 <- pbmc[["apca"]]@cell.embeddings

svd_1 <- irlba::irlba(embedding_1, nv = 3)
embedding_1 <- embedding_1/svd_1$d[1]*nrow(embedding_1)
svd_2 <- irlba::irlba(embedding_2, nv = 3)
embedding_2 <- embedding_2/svd_2$d[1]*nrow(embedding_2)

set.seed(10)
both_embedding <- Seurat::RunUMAP(cbind(embedding_1, embedding_2))
both_embedding <- both_embedding@cell.embeddings
pbmc[["both"]] <- Seurat::CreateDimReducObject(embedding = both_embedding, key = "both_", assay = "SCT")

###########################


title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

# plot RNA embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(pbmc, reduction = main_vec[i],
                           group.by = "celltype.l2", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("CITE-seq PBMC (RNA)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_citeseq_pbmc_dcca_rna_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot ATAC embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(pbmc, reduction = paste0(main_vec[i], "2"),
                           group.by = "celltype.l2", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("CITE-seq PBMC (Protein)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_citeseq_pbmc_dcca_protein_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot the combined view
plot1 <- Seurat::DimPlot(pbmc, reduction = "combined",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("CITE-seq PBMC\nBoth Common")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_citeseq_pbmc_dcca_both_common_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# plot the everything combined for both
plot1 <- Seurat::DimPlot(pbmc, reduction = "both",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("CITE-seq PBMC\nBoth Everything")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_citeseq_pbmc_dcca_both_everything_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


#########################

# plot the Seurat embeddings
plot1 <- Seurat::DimPlot(pbmc, reduction = "rna.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("CITE-seq PBMC (RNA)\nSeurat UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_citeseq_pbmc_dcca_seurat_rna.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "adt.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("CITE-seq PBMC (Protein)\nSeurat UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_citeseq_pbmc_dcca_seurat_protein.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "wnn.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("CITE-seq PBMC (WNN)\nSeurat UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_citeseq_pbmc_dcca_seurat_wnn.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


#######################

png("../../../../out/figures/Writeup14d/Writeup14d_citeseq_pbmc_dcca_summary.png",
    height = 1500, width = 1500, units = "px", res = 300)
par(mar = c(5,5,5,5))
multiomicCCA::plot_summary(dcca_res,
                           main = "CITE-seq PBMC")
graphics.off()

png("../../../../out/figures/Writeup14d/Writeup14d_citeseq_pbmc_dcca_scores.png",
    height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap.dcca(dcca_res,
                                       membership_vec = as.factor(pbmc@meta.data$celltype.l2),
                                       log_scale = T, scaling_power = 4)
graphics.off()
