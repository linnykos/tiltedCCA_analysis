rm(list=ls())

library(Seurat); library(Signac); library(slingshot)

load("../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")
greenleaf_tmp <- greenleaf
load("../../../out/main/10x_greenleaf_preprocessed.RData")
greenleaf[["common_tcca"]] <- greenleaf_tmp[["common_tcca"]]

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()


Seurat::DefaultAssay(greenleaf) <- "SCT"
set.seed(10)
greenleaf <- Seurat::FindNeighbors(greenleaf, dims = 1:50, reduction = "lsi")
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf, resolution = 2)

Seurat::DefaultAssay(greenleaf) <- "ATAC"
set.seed(10)
greenleaf <- Seurat::FindNeighbors(greenleaf, dims = 2:50, reduction = "lsi")
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf, resolution = 2)

Seurat::DefaultAssay(greenleaf) <- "ATAC"
set.seed(10)
greenleaf <- Seurat::FindNeighbors(greenleaf, dims = 2:50, reduction = "lsi")
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf, resolution = 2)
