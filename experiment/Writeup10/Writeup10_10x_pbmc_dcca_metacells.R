rm(list=ls()); set.seed(10); # gcinfo(TRUE)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../out/Writeup10_10x_pbmc_preprocess4.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

#############

mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ATAC"]]@scale.data)

set.seed(10)
rank_1 <- 17; rank_2 <- 13; K <- min(c(rank_1, rank_2))
meta_clustering <- stats::kmeans(pbmc[["umap.rna"]]@cell.embeddings, centers = 50*K)$cluster
table(table(meta_clustering))

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
rm(list = "pbmc"); gc(T)

gcinfo(F)
set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, rank_1 = rank_1, rank_2 = rank_2,
                                      meta_clustering = meta_clustering,
                                      apply_shrinkage = F, verbose = T) # takes around 8 minutes
rm(list = c("mat_1", "mat_2"))

source_code <- readLines("../multiomicCCA_analysis/experiment/Writeup10/Writeup10_10x_pbmc_dcca_metacells.R")
save.image("../../out/Writeup10_10x_pbmc_dcca_metacells2.RData")





