rm(list=ls())
library(Seurat); library(Signac)

source("greenleaf_colorPalette.R")
load("../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")
# load("../../../out/Writeup14n/10x_greenleaf_tcca_RNA-ATAC.RData")

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

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "celltype", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_RNA-ATA_cleaned-smaller.png"),
                plot1, device = "png", width = 2, height = 2, units = "in",
                dpi = 500)


#####################

zz <- greenleaf[["common_tcca"]]@cell.embeddings
zz <- zz[-intersect(which(zz[,1] <= 0), which(zz[,2]>0)),]
zz <- zz[-which(zz[,2] >= 10),]
zz <- zz[-which(zz[,2] <= -9),]
nrow(zz)

