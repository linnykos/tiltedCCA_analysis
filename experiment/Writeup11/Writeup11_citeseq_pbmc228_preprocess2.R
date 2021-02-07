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

#######################

pbmc[["spca"]] <- NULL
pbmc[["umap"]] <- NULL
pbmc[["aumap"]] <- NULL
pbmc[["wknn"]] <- NULL
pbmc[["wsnn"]] <- NULL
pbmc[["wnn.umap"]] <- NULL

pbmc[["SCT"]]@scale.data <- tcrossprod(pbmc[["pca"]]@feature.loadings, pbmc[["pca"]]@cell.embeddings)
pbmc[["ADT"]]@scale.data <- tcrossprod(pbmc[["apca"]]@feature.loadings, pbmc[["apca"]]@cell.embeddings)
dim(pbmc[["SCT"]]@scale.data)
dim(pbmc[["ADT"]]@scale.data)

set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'pca', dims = 1:40, assay = 'SCT',
                        reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'apca', dims = 1:50, assay = 'ADT',
                        reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

# the following takes a hot minute with 160K cells
pbmc <- Seurat::FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "apca"),
                                        dims.list = list(1:40, 1:50), modality.weight.name = "RNA.weight")
head(pbmc@meta.data)
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

save.image("../../../../out/Writeup11/Writeup11_citeseq_pbmc228_preprocessed2.RData")


