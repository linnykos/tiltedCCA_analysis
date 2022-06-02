rm(list=ls())

library(Seurat); library(Signac); library(Slingshot)

load("../../../out/main/10x_greenleaf_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(greenleaf) <- "ATAC"
set.seed(10)
greenleaf <- Seurat::FindNeighbors(greenleaf, dims = 1:50, reduction = "lsi")
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf, resolution = 2)
