rm(list=ls())
library(Seurat); library(Signac)

source("greenleaf_colorPalette.R")
load("../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = FALSE,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-ATAC_umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::FeaturePlot(greenleaf, reduction = "common_tcca",
                             features = "pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC): T-CCA's Common\nSlingshot's pseudotime via ATAC"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-ATAC_umap_common-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(greenleaf, reduction = "distinct1_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = FALSE,
                         cols = col_palette)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nRNA distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-ATAC_umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(greenleaf, reduction = "distinct2_tcca",
                         group.by = "celltype", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster = FALSE,
                         cols = col_palette)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nATAC distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-ATAC_umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

############

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap",
                         group.by = "celltype", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_rna-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap.atac",
                         group.by = "celltype", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_atac-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "celltype", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-ATA_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap.wnn",
                         group.by = "celltype", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_wnn-umap_RNA-ATA_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "consensusUMAP",
                         group.by = "celltype", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_consensusPCA-umap_RNA-ATA_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "celltype", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-ATA_cleaned-smaller.png"),
                plot1, device = "png", width = 2, height = 2, units = "in",
                dpi = 500)

#################################################

keep_vec <- rep(1, ncol(greenleaf))
keep_vec[greenleaf$celltype %in% c("EC/Peric.", "IN1", "IN2", "IN3", "mGPC/OPC", "SP")] <- 0
greenleaf$keep <- keep_vec
greenleaf <- subset(greenleaf, keep == 1)

greenleaf[["umap.atac"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = greenleaf[["umap"]]@cell.embeddings,
  target_mat = greenleaf[["umap.atac"]]@cell.embeddings
)

keep_vec <- rep(1, ncol(greenleaf))
keep_vec[greenleaf[["umap"]]@cell.embeddings[,2]>=5] <- 0
keep_vec[greenleaf[["umap"]]@cell.embeddings[,2]<=-11] <- 0
keep_vec[greenleaf[["umap.atac"]]@cell.embeddings[,2]>=5.5] <- 0
keep_vec[greenleaf[["umap.atac"]]@cell.embeddings[,2]<=-11] <- 0
greenleaf$keep <- keep_vec
greenleaf <- subset(greenleaf, keep == 1)

greenleaf[["umap.atac"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = greenleaf[["umap"]]@cell.embeddings,
  target_mat = greenleaf[["umap.atac"]]@cell.embeddings
)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap.atac",
                         group.by = "celltype", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("ATAC UMAP 2") + ggplot2::xlab("ATAC UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_atac-umap_slides.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

plot1 <- Seurat::DimPlot(greenleaf, reduction = "umap",
                         group.by = "celltype", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("RNA UMAP 2") + ggplot2::xlab("RNA UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_rna-umap_slides.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")


greenleaf[["common_tcca"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = greenleaf[["umap"]]@cell.embeddings,
  target_mat = greenleaf[["common_tcca"]]@cell.embeddings
)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "celltype", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Common UMAP 2") + ggplot2::xlab("Common UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-ATAC_common-umap_slides.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")


