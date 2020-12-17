rm(list=ls())

# from https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html

library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)

date_of_run <- Sys.time()
session_info <- sessionInfo()

SeuratData::InstallData("bmcite")
bm <- SeuratData::LoadData(ds = "bmcite")

DefaultAssay(bm) <- 'RNA'
bm <- Seurat::NormalizeData(bm) %>% Seurat::FindVariableFeatures() %>% Seurat::ScaleData() %>% Seurat::RunPCA()

DefaultAssay(bm) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
Seurat::VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- Seurat::NormalizeData(bm, normalization.method = 'CLR', margin = 2) %>% 
  Seurat::ScaleData() %>% Seurat::RunPCA(reduction.name = 'apca')

##########################

bm <- Seurat::RunUMAP(bm, reduction = 'pca', dims = 1:30, assay = 'RNA', 
              reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
bm <- Seurat::RunUMAP(bm, reduction = 'apca', dims = 1:18, assay = 'ADT', 
              reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

names(bm)

pca_log_rna <- bm[["pca"]]@cell.embeddings[,1:20]
pca_clr_adt <- bm[["apca"]]@cell.embeddings[,1:20]
umap_log_rna <- bm[["rna.umap"]]@cell.embeddings
umap_clr_adt <- bm[["adt.umap"]]@cell.embeddings

save(pca_log_rna, pca_clr_adt, umap_log_rna, umap_clr_adt, 
     file = "../../out/Writeup9_citeseq_bonemarrow_dimred_onlyembedding.RData")

source_code <- readLines("Writeup9_citeseq_bonemarrow_dimred.R")
save.image(file = "../../out/Writeup9_citeseq_bonemarrow_dimred_all.RData")



