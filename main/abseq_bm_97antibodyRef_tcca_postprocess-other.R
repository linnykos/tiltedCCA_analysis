rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/abseq_bm97Ref_preprocessed.RData")
source("bm_97antibodyRef_colorPalette-simpler.R")

bm[["adt.umap"]]@cell.embeddings <- tiltedCCA:::.rotate_matrix(
  source_mat = bm[["rna.umap"]]@cell.embeddings,
  target_mat = bm[["adt.umap"]]@cell.embeddings
)
plot1 <- Seurat::DimPlot(bm, reduction = "adt.umap",
                         group.by = "ct", 
                         cols = col_palette, pt.size = .25)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("AB UMAP 2") + ggplot2::xlab("AB UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_adt-umap_slides-simpler.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

plot1 <- Seurat::DimPlot(bm, reduction = "rna.umap",
                         group.by = "ct", 
                         cols = col_palette, pt.size = .25)
plot1 <- plot1 + Seurat::NoLegend() 
plot1 <- plot1 + ggplot2::ggtitle("") + ggplot2::ylab("RNA UMAP 2") + ggplot2::xlab("RNA UMAP 1")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_rna-umap_slides-simpler.png"),
                plot1, device = "png", width = 4, height = 4, units = "in")

