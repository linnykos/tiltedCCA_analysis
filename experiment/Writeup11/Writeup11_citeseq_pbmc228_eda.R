rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_pbmc228_preprocessed.RData")

###########################
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)

#3########################

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc_rna_donor_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'rna.umap', group.by = 'donor', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA UMAP (Donor, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc_protein_donor_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'adt.umap', group.by = 'donor', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Protein UMAP (Donor, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc_rna_time_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'rna.umap', group.by = 'time', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA UMAP (Time, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc_time_donor_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'adt.umap', group.by = 'time', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Protein UMAP (Time, Vanilla)")
graphics.off()