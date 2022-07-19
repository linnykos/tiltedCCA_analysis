rm(list=ls())
library(Seurat)

load("../../../out/main/10x_reik_tiltedcca.RData")
source("reik_colorPalette.R")

plot1 <- Seurat::DimPlot(reik, reduction = "common_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_tcca-umap_common.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "common_tcca",
                         group.by = "orig.ident", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_tcca-umap_common_orig-ident.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "distinct1_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nRNA Distinct subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_tcca-umap_distinct1.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")

plot1 <- Seurat::DimPlot(reik, reduction = "distinct2_tcca",
                        group.by = "celltype", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Mouse embryo (10x, RNA+ATAC)\nATAC Distinct subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_tcca-umap_distinct2.png"),
                plot1, device = "png", width = 10, height = 5, units = "in")


