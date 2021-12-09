rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")
source("../../simulation/scai.R")
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat_1 <- t(as.matrix(bm[["RNA"]]@data))
mat_1 <- mat_1[,which(colnames(mat_1) %in% rownames(bm[["RNA"]]@scale.data))]
mat_2 <- t(as.matrix(bm[["ADT"]]@data))
mat_2 <- mat_2[,which(colnames(mat_2) %in% rownames(bm[["ADT"]]@scale.data))]
scai_res <- scai(mat_1, mat_2, r = min(ncol(mat_1), ncol(mat_2)), gamma = 0)
save.image(file = "../../../../out/Writeup14g/Writeup14g_citeseq_bm25_scai.RData")


