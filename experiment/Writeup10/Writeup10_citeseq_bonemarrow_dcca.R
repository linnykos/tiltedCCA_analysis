rm(list=ls())

library(Seurat)
load("../../out/Writeup9_citeseq_bonemarrow_dimred_all.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- sessionInfo()

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)

library(multiomicCCA)

set.seed(10)
K <- 8
dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K, verbose = T) 

rm(list = c("mat_1", "mat_2", "bm", "pca_clr_adt", "pca_log_rna", "umap_clr_adt", "umap_log_rna"))

source_code <- readLines("../multiomicCCA_analysis/experiment/Writeup10/Writeup10_citeseq_bonemarrow_dcca.R")
save.image("../../out/Writeup10_10x_pbmc_dcca.RData")
