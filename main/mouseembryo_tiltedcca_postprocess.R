rm(list=ls())
library(Seurat); library(Signac)

source("mouseembryo_colorPalette.R")
load("../../../out/main/10x_mouseembryo_tiltedcca.RData")

plot1 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE,
                         cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_tcca-umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::FeaturePlot(mbrain, reduction = "common_tcca",
                             features = "pseudotime")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC): T-CCA's Common\nSlingshot's pseudotime via ATAC"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_tcca-umap_common-pseudotime.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(mbrain, reduction = "distinct1_tcca",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE,
                         cols = col_palette)
plot2 <- plot2 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nRNA distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_tcca-umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(mbrain, reduction = "distinct2_tcca",
                         group.by = "label_Savercat", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE,
                         cols = col_palette)
plot3 <- plot3 + ggplot2::ggtitle(paste0("Mouse Embryo E18 (10x, RNA+ATAC)\nATAC distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_tcca-umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

#########################

plot1 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         group.by = "label_Savercat", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_tcca-umap_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)


plot1 <- Seurat::DimPlot(mbrain, reduction = "common_tcca",
                         group.by = "label_Savercat", 
                         cols = col_palette)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_mouseembryo_tcca-umap_cleaned-smaller.png"),
                plot1, device = "png", width = 2, height = 2, units = "in",
                dpi = 500)


