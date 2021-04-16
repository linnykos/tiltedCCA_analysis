rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)

set.seed(10)
rank_1 <- 30; rank_2 <- 18
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, rank_1 = rank_1, rank_2 = rank_2, 
                                      apply_shrinkage = F, verbose = T) 
dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, rank_c = min(rank_1, rank_2), verbose = T)
svd_list <- multiomicCCA::extract_svd_embedding(dcca_decomp)
