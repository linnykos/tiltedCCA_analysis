rm(list=ls())
library(Seurat)

load("../../../../out/Writeup14n/abseq_bm97Ref_tcca.RData")
source("../../main/bm_97antibodyRef_colorPalette.R")

plot1 <-Seurat::DimPlot(bm, reduction = "common_tcca",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14n/abseq_bm97Ref_tcca-umap_common.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "distinct1_tcca",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nRNA Distinct subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14n/abseq_bm97Ref_tcca-umap_distinct1.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "distinct2_tcca",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nADT Distinct subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14n/abseq_bm97Ref_tcca-umap_distinct2.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")
