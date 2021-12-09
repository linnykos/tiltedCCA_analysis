rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

bm2 <- bm
Seurat::DefaultAssay(bm2) <- "RNA"
bm2[["ADT"]] <- NULL
bm2[["spca"]] <- NULL; bm2[["adt.umap"]] <- NULL; bm2[["wnn.umap"]] <- NULL
mat <- bm[["ADT"]]@counts
protein_idx <- which(!rownames(mat) %in% c("CD8a", "CD4"))