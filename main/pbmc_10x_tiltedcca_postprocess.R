rm(list=ls())
load("../../../out/main/10x_pbmc_tiltedcca.RData")

library(Seurat); library(Signac)

plot1 <- Seurat::DimPlot(pbmc, reduction = "common_tcca",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_pbmc_tcca-umap_common.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot2 <- Seurat::DimPlot(pbmc, reduction = "distinct1_tcca",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("PBMC (10x, RNA+ATAC)\nRNA distinct subspace"))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_pbmc_tcca-umap_distinct1.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(pbmc, reduction = "distinct2_tcca",
                         group.by = "predicted.id", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("PBMC (10x, RNA+ATAC)\nATAC distinct subspace"))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_pbmc_tcca-umap_distinct2.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")


