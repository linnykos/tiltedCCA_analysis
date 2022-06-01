rm(list=ls())
load("../../../../out/main/10x_greenleaf_preprocessed.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(greenleaf)

Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
Seurat::DefaultAssay(greenleaf) <- "ATAC"
mat_2 <- Matrix::t(greenleaf[["ATAC"]]@data[Seurat::VariableFeatures(object = greenleaf),])

mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}
