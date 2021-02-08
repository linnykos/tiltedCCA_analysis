rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_pbmc228_preprocessed2.RData")

###########################
library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)

#########
# one-time check:
zz <- Seurat::RunPCA(pbmc[["SCT"]]@scale.data)
# looks roughly the same
zz@cell.embeddings[1:10,1:5]
pbmc[["pca"]]@cell.embeddings[1:10,1:5]
rm(list = "zz")

#3########################

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc224_rna_custom_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'rna.umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA UMAP (224, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc224_protein_custom_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'adt.umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Protein UMAP (224, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc224_wnn_custom_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("WNN UMAP (224, Seurat)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc224_seurat_weight.png", height = 1500, width = 2500, units = "px", res = 300)
Seurat::VlnPlot(pbmc, features = "RNA.weight", group.by = 'celltype.l2', sort = TRUE, pt.size = 0.1) +
  NoLegend() + ggplot2::ggtitle("Seurat WNN's RNA weight (224)")
graphics.off()

#############################

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc224_rna_donor_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'rna.umap', group.by = 'donor', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA UMAP (Donor, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc224_protein_donor_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'adt.umap', group.by = 'donor', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Protein UMAP (Donor, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc224_wnn_donor_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'wnn.umap', group.by = 'donor', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("WNN UMAP (Donor, Seurat)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc224_rna_time_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'rna.umap', group.by = 'time', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA UMAP (Time, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc224_protein_time_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'adt.umap', group.by = 'time', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Protein UMAP (Time, Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc224_wnn_time_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'wnn.umap', group.by = 'time', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("WNN UMAP (Time, Seurat)")
graphics.off()
