rm(list=ls())
library(Seurat)

load("../../../out/main/abseq_bm97Ref_tcca.RData")
source("bm_97antibodyRef_colorPalette.R")

plot1 <-Seurat::DimPlot(bm, reduction = "common_tcca",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_tcca-umap_common.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "distinct1_tcca",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nRNA Distinct subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_tcca-umap_distinct1.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

plot1 <-Seurat::DimPlot(bm, reduction = "distinct2_tcca",
                        group.by = "ct", label = TRUE,
                        repel = TRUE, label.size = 2.5,
                        cols = col_palette)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human BM (Abseq, RNA+ADT)\nADT Distinct subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_tcca-umap_distinct2.png"),
                plot1, device = "png", width = 11, height = 5, units = "in")

#############################


bm[["common_tcca"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = bm[["rna.umap"]]@cell.embeddings,
  target_mat = bm[["common_tcca"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(bm, reduction = "common_tcca",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 1.5)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Common UMAP 2") + ggplot2::xlab("Common UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_tcca-umap_common_slides.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")


bm[["distinct2_tcca"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = bm[["rna.umap"]]@cell.embeddings,
  target_mat = bm[["distinct2_tcca"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(bm, reduction = "distinct2_tcca",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 0.05)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Distinct UMAP 2") + ggplot2::xlab("Distinct UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_tcca-umap_distinct2_slides.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

bm[["common_tcca"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = bm[["rna.umap"]]@cell.embeddings,
  target_mat = bm[["common_tcca"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(bm, reduction = "common_tcca",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 0.05)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("Common UMAP 2") + ggplot2::xlab("Common UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_tcca-umap_common_slides-small.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

bm[["adt.umap"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = bm[["rna.umap"]]@cell.embeddings,
  target_mat = bm[["adt.umap"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 0.05)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("AB UMAP 2") + ggplot2::xlab("AB UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_adt-umap_slides-small.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                         group.by = "ct", 
                         cols = col_palette, pt.size = 0.05)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("RNA UMAP 2") + ggplot2::xlab("RNA UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_rna-umap_slides-small.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")



