rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
source("../Writeup14/Writeup14_peakcalling_function.R")
date_of_run <- Sys.time(); session_info <- sessionInfo()

mat_1 <- t(mbrain[["SCT"]]@scale.data)
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- Matrix::t(mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),])

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- mbrain@meta.data

#############


cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Hindbrain glycinergic", "Midbrain glutamatergic",
                                                 "Forebrain glutamatergic", "Radial glia", "Forebrain GABAergic", 
                                                 "Neuroblast", "Cajal-Retzius", "Mixed region GABAergic", 
                                                 "Glioblast", "Cortical or hippocampal glutamatergic"))
set.seed(10)
rank_1 <- 30; rank_2 <- 50
mat_1 <- mat_1[cell_idx,]; mat_2 <- mat_2[cell_idx,]
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = T,
                                      meta_clustering = NA, num_neigh = 15, 
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 

##################

membership_vec <- as.factor(metadata$label_Savercat[cell_idx])

set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = 15, membership_vec = membership_vec,
                                         data_1 = T, data_2 = F,
                                         radius_quantile = 0.5,
                                         bool_matrix = T, verbose = T)

set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = 15, membership_vec = membership_vec,
                                          data_1 = F, data_2 = T,
                                          radius_quantile = 0.5,
                                          bool_matrix = T, verbose = T)

g <- atac_frnn$c_g
zz <- sapply(1:nrow(g), function(x){length(multiomicCCA:::.nonzero_col(g, x, bool_value = F))})
quantile(zz)
g2 <- Matrix::t(g)
zz <- sapply(1:nrow(g2), function(x){length(multiomicCCA:::.nonzero_col(g2, x, bool_value = F))})
quantile(zz)
g <- multiomicCCA:::.symmetrize_sparse(g, set_ones = F)
zz <- sapply(1:nrow(g), function(x){length(multiomicCCA:::.nonzero_col(g, x, bool_value = F))})
quantile(zz)
