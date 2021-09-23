rm(list=ls())
library(Seurat); library(Signac)
load("../../../../out/Writeup14d/Writeup14d_10x_pbmc_dcca.RData")

pbmc[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "common_", assay = "RNA")
pbmc[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "distinct_", assay = "RNA")
pbmc[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "everything_", assay = "RNA")
pbmc[["clap"]] <- Seurat::CreateDimReducObject(embedding = c_eig, key = "clap_", assay = "RNA")
pbmc[["dlap"]] <- Seurat::CreateDimReducObject(embedding = d_eig, key = "dlap_", assay = "RNA")
pbmc[["elap"]] <- Seurat::CreateDimReducObject(embedding = e_eig, key = "elap_", assay = "RNA")
rm(list = c("c_eig", "d_eig", "e_eig", "rna_embeddings", "rna_frnn")); gc(T)

pbmc[["common2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[1]], key = "common2_", assay = "RNA")
pbmc[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[2]], key = "distinct2_", assay = "RNA")
pbmc[["everything2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[3]], key = "everything2_", assay = "RNA")
pbmc[["clap2"]] <- Seurat::CreateDimReducObject(embedding = c_eig2, key = "clap2_", assay = "RNA")
pbmc[["dlap2"]] <- Seurat::CreateDimReducObject(embedding = d_eig2, key = "dlap2_", assay = "RNA")
pbmc[["elap2"]] <- Seurat::CreateDimReducObject(embedding = e_eig2, key = "elap2_", assay = "RNA")
rm(list = c("c_eig2", "d_eig2", "e_eig2", "atac_embeddings", "atac_frnn")); gc(T)

pbmc[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, key = "combined_", assay = "RNA")
rm(list = c("combined_common_umap")); gc(T)

embedding_1 <- pbmc[["pca"]]@cell.embeddings
embedding_2 <- pbmc[["lsi"]]@cell.embeddings
set.seed(10)
both_embedding <- Seurat::RunUMAP(cbind(embedding_1, embedding_2))
both_embedding <- both_embedding@cell.embeddings
pbmc[["both"]] <- Seurat::CreateDimReducObject(embedding = both_embedding, key = "both_", assay = "RNA")

###########################

# done with all the calculations (for now). now plot

title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

# plot RNA embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(pbmc, reduction = main_vec[i],
                           group.by = "predicted.id", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("10x PBMC (RNA)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_10x_pbmc_dcca_rna_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot ATAC embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(pbmc, reduction = paste0(main_vec[i], "2"),
                           group.by = "predicted.id", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("10x PBMC (ATAC)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_10x_pbmc_dcca_atac_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# plot the combined view
plot1 <- Seurat::DimPlot(pbmc, reduction = "combined",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("10x PBMC\nBoth Common")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_10x_pbmc_dcca_both_common_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

# plot the everything combined for both
plot1 <- Seurat::DimPlot(pbmc, reduction = "both",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("10x PBMC\nBoth Everything")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_10x_pbmc_dcca_both_everything_umap.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


#########################

# plot the Seurat embeddings
plot1 <- Seurat::DimPlot(pbmc, reduction = "umap.rna",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("10x PBMC (RNA)\nSeurat UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_10x_pbmc_dcca_seurat_rna.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "umap.atac",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("10x PBMC (ATAC)\nSeurat UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_10x_pbmc_dcca_seurat_atac.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(pbmc, reduction = "wnn.umap",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("10x PBMC (WNN)\nSeurat UMAP"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_10x_pbmc_dcca_seurat_wnn.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

