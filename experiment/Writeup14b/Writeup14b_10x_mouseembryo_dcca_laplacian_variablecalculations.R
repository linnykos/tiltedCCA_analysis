rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

mat_1 <- t(mbrain[["SCT"]]@scale.data)
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- Matrix::t(mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),])

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- mbrain@meta.data

rm(list = c("mbrain"))
gc()

#########

cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Hindbrain glycinergic", "Midbrain glutamatergic",
                                                 "Forebrain glutamatergic", "Radial glia", "Forebrain GABAergic", 
                                                 "Neuroblast", "Cajal-Retzius", "Mixed region GABAergic", 
                                                 "Glioblast", "Cortical or hippocampal glutamatergic"))
set.seed(10)
rank_1 <- 30; rank_2 <- 31
mat_1 <- mat_1[cell_idx,]; mat_2 <- mat_2[cell_idx,]
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = T,
                                      meta_clustering = NA, num_neigh = 15, 
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
mat_2_denoised <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2
membership_vec <- as.factor(metadata$label_Savercat[cell_idx])

rm(list = c("mat_2"))
gc()

save(date_of_run, session_info, dcca_res, membership_vec,
     file = "../../../../out/Writeup14b/Writeup14b_10x_mouseembryo_dcca_laplacian_variablecalculations.RData")

###########

# do the RNA ones 
set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = 15, membership_vec = membership_vec,
                                         data_1 = T, data_2 = F,
                                         bool_matrix = T, include_diag = F, verbose = T)

#compute all the degree vectors
k_max <- 200
c_eig <- multiomicCCA::compute_laplacian(rna_frnn$c_g, k_max = k_max, rowname_vec = rownames(mat_1), 
                                         colname_vec = paste0("clap_", 1:k_max))
d_eig <- multiomicCCA::compute_laplacian(rna_frnn$d_g, k_max = k_max, rowname_vec = rownames(mat_1), 
                                         colname_vec = paste0("dlap_", 1:k_max))
e_eig <- multiomicCCA::compute_laplacian(rna_frnn$e_g, k_max = k_max, rowname_vec = rownames(mat_1), 
                                         colname_vec = paste0("elap_", 1:k_max))

p1 <- ncol(mat_1)
gene_smoothed <- lapply(1:p1, function(j){
  if(j %% floor(p1/10) == 0) cat('*')
  
  c_res <- compute_smooth_signal(mat_1_denoised[,j], c_eig)
  d_res <- compute_smooth_signal(mat_1_denoised[,j], d_eig)
  e_res <- compute_smooth_signal(mat_1_denoised[,j], e_eig)
  
  list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
       d_variance = d_res$variance, d_r2 = d_res$r_squared,
       e_variance = e_res$variance, e_r2 = e_res$r_squared)
})

save(date_of_run, session_info, dcca_res, membership_vec, 
     rna_frnn, c_eig, d_eig, e_eig, gene_smoothed,
     file = "../../../../out/Writeup14b/Writeup14b_10x_mouseembryo_dcca_laplacian_variablecalculations.RData")

###############################3

set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = 15, membership_vec = membership_vec,
                                          data_1 = F, data_2 = T,
                                          bool_matrix = T, include_diag = F, verbose = T)

k_max <- 200
c_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$c_g, k_max = k_max, rowname_vec = rownames(mat_1), 
                                          colname_vec = paste0("clap_", 1:k_max))
d_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$d_g, k_max = k_max, rowname_vec = rownames(mat_1), 
                                          colname_vec = paste0("dlap_", 1:k_max))
e_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$e_g, k_max = k_max, rowname_vec = rownames(mat_1), 
                                          colname_vec = paste0("elap_", 1:k_max))

p2 <- ncol(mat_2_denoised)
atac_smoothed <- lapply(1:p2, function(j){
  print(p2)
  
  c_res <- compute_smooth_signal(mat_2_denoised[,j], c_eig)
  d_res <- compute_smooth_signal(mat_2_denoised[,j], d_eig)
  e_res <- compute_smooth_signal(mat_2_denoised[,j], e_eig)
  
  list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
       d_variance = d_res$variance, d_r2 = d_res$r_squared,
       e_variance = e_res$variance, e_r2 = e_res$r_squared)
})

save(date_of_run, session_info, dcca_res, membership_vec, 
     rna_frnn, c_eig, d_eig, e_eig, gene_smoothed,
     atac_frnn, c_eig2, d_eig2, e_eig2, atac_smoothed,
     file = "../../../../out/Writeup14b/Writeup14b_10x_mouseembryo_dcca_laplacian_variablecalculations.RData")


