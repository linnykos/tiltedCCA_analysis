rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

###########################
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)

png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_rna_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'rna.umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA UMAP (25, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_protein_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'adt.umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Protein UMAP (25, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_wnn_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(bm, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("WNN UMAP (25, Seurat)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_bm25_seurat_weight.png", height = 1500, width = 2500, units = "px", res = 300)
Seurat::VlnPlot(bm, features = "RNA.weight", group.by = 'celltype.l2', sort = TRUE, pt.size = 0.1) +
  NoLegend() + ggplot2::ggtitle("Seurat WNN's RNA weight (25)")
graphics.off()