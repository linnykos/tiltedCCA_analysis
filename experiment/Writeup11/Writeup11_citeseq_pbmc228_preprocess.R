# rm(list=ls())
# library(Seurat)
# library(SeuratDisk)
# library(SeuratData)
# 
# pbmc <- SeuratDisk::LoadH5Seurat("../../../../data/CITE_PBMC_RNA-228Protein/multi.h5seurat")
# pbmc
# dim(pbmc[["SCT"]])

##############
rm(list=ls())
load("../../../../data/CITE_PBMC_RNA-228Protein/pbmc.RData")
names(pbmc)

library(Seurat)
library(SeuratDisk)
library(SeuratData)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(dplyr)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)


png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc_rna_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("RNA UMAP (Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc_protein_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'aumap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("Protein UMAP (Vanilla)")
graphics.off()

png("../../../../out/figures/Writeup11/Writeup11_citeseq_pbmc_wnn_umap.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pbmc, reduction = 'wnn.umap', group.by = 'celltype.l2', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("WNN UMAP (Seurat)")
graphics.off()

##############

pbmc[["wknn"]] <- NULL
pbmc[["wsnn"]] <- NULL
pbmc[["apca"]] <- NULL
pbmc[["aumap"]] <- NULL
pbmc[["pca"]] <- NULL
pbmc[["spca"]] <- NULL
pbmc[["umap"]] <- NULL
pbmc[["wnn.umap"]] <- NULL

# preprocess data
Seurat::DefaultAssay(pbmc) <- "SCT"
dim(pbmc)
pbmc <- Seurat::ScaleData(pbmc, split.by = list(pbmc@meta.data$donor, pbmc@meta.data$time))
pbmc <- Seurat::RunPCA(pbmc)
dim(pbmc[["SCT"]]@scale.data)

Seurat::DefaultAssay(pbmc) <- "ADT"
dim(pbmc)
Seurat::VariableFeatures(pbmc) <- rownames(pbmc[["ADT"]])
pbmc <- Seurat::ScaleData(pbmc, split.by = list(pbmc@meta.data$donor, pbmc@meta.data$time))
pbmc <- Seurat::RunPCA(pbmc, reduction.name = 'apca')
dim(pbmc[["ADT"]]@scale.data)

head(pbmc@meta.data)
table(pbmc@meta.data[,"celltype.l2"])
names(pbmc)

# the following takes a hot minute with 160K cells
pbmc <- Seurat::FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "apca"),
                                        dims.list = list(1:40, 1:50), modality.weight.name = "RNA.weight")
head(pbmc@meta.data)
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

pbmc[["pca"]]@stdev
abs(diff(pbmc[["pca"]]@stdev)/pbmc[["pca"]]@stdev[-1])
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'pca', dims = 1:40, assay = 'SCT',
                        reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'apca', dims = 1:50, assay = 'ADT',
                        reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')


save.image("../../../../out/Writeup11/Writeup11_citeseq_pbmc228_preprocessed.RData")

