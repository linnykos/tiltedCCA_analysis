rm(list=ls())
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")
source("pbmc_224antibody_colorPalette.R")

library(Seurat); library(Signac)

plot1 <- Seurat::DimPlot(pbmc, reduction = "common_tcca",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(pbmc, reduction = "distinct1_tcca",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nRNA distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(pbmc, reduction = "distinct2_tcca",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nADT distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

###############################

plot1 <- Seurat::DimPlot(pbmc, reduction = "common_tcca",
                         group.by = "celltype.l2",
                         cols = col_palette,
                         raster = FALSE)
plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
plot1 <- plot1 + ggplot2::ggtitle("")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_common_cleaned.png"),
                plot1, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot2 <- Seurat::DimPlot(pbmc, reduction = "distinct1_tcca",
                         group.by = "celltype.l2", 
                         cols = col_palette,
                         raster = FALSE)
plot2 <- plot2 + Seurat::NoLegend() + Seurat::NoAxes()
plot2 <- plot2 + ggplot2::ggtitle("")
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_distinct1_cleaned.png"),
                plot2, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot3 <- Seurat::DimPlot(pbmc, reduction = "distinct2_tcca",
                         group.by = "celltype.l2",
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + Seurat::NoLegend() + Seurat::NoAxes()
plot3 <- plot3 + ggplot2::ggtitle("")
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_distinct2_cleaned.png"),
                plot3, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)



