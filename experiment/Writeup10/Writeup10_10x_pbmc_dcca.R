rm(list=ls()); set.seed(10); gcinfo(TRUE)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../out/Writeup10_10x_pbmc_preprocess4.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

#############

mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ATAC"]]@scale.data)

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
rm(list = "pbmc"); gc(T)

set.seed(10)
K <- 7
dcca_res <- dcca_factor(mat_1, mat_2, rank_1 = K, rank_2 = K) # takes around 20 min

source_code <- readLines("experiment/Writeup10_10x_pbmc_dcca.R")
save.image("../../out/Writeup10_10x_pbmc_dcca.RData")
