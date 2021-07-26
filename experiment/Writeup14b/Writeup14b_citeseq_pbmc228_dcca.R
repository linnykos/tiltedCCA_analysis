rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_pbmc228_preprocessed2.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ADT"]]@scale.data)
dim(mat_1); dim(mat_2)
metadata <- pbmc@meta.data
rm(list = "pbmc"); gc(T)
membership_vec <- as.factor(metadata$celltype.l2)

set.seed(10)
idx <- multiomicCCA::construct_celltype_subsample(membership_vec, min_subsample_cell = 3000)
membership_vec <- membership_vec[idx]
mat_1 <- mat_1[idx,]; mat_2 <- mat_2[idx,]

set.seed(10)
rank_1 <- 40; rank_2 <- 50; nn <- 30
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      meta_clustering = NA, num_neigh = nn, cell_max = 15000,
                                      apply_shrinkage = F, fix_distinct_perc = F, 
                                      verbose = T) 
cell_name <- rownames(mat_1)

save(date_of_run, session_info, dcca_res,
     file = "../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca1.RData")

rm(list = c("mat_1", "mat_2")); gc(T)
##################

set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                         data_1 = T, data_2 = F,
                                         radius_quantile = 0.5, symmetrize = F, 
                                         bool_matrix = T, verbose = T)

save(date_of_run, session_info, rna_frnn, 
     file = "../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca2.RData")

#compute all rna basis vectors
k_max <- 50
c_eig <- multiomicCCA::compute_laplacian(rna_frnn$c_g, k_max = k_max, rowname_vec = cell_name, 
                                         colname_vec = paste0("clap_", 1:k_max), verbose = F)
d_eig <- multiomicCCA::compute_laplacian(rna_frnn$d_g, k_max = k_max, rowname_vec = cell_name, 
                                         colname_vec = paste0("dlap_", 1:k_max), verbose = F)
e_eig <- multiomicCCA::compute_laplacian(rna_frnn$e_g, k_max = k_max, rowname_vec = cell_name, 
                                         colname_vec = paste0("elap_", 1:k_max), verbose = F)

print("Finished eigenbasis")

set.seed(10)
rna_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = T, data_2 = F,
                                                 c_g = rna_frnn$c_g, d_g = rna_frnn$d_g, 
                                                 only_embedding = T, verbose = T, 
                                                 sampling_type = "adaptive_gaussian",
                                                 keep_nn = T)

save(date_of_run, session_info, rna_frnn, c_eig, d_eig, e_eig, rna_embeddings, 
     file = "../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca2.RData")

rm(list = "rna_embeddings"); gc(T)

################

set.seed(10)
protein_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                             data_1 = F, data_2 = T,
                                             radius_quantile = 0.5,
                                             bool_matrix = T, symmetrize = F, verbose = T)

#compute all the degree vectors
k_max <- 50
c_eig2 <- multiomicCCA::compute_laplacian(protein_frnn$c_g, k_max = k_max, rowname_vec = cell_name, 
                                          colname_vec = paste0("clap_", 1:k_max), verbose = F)
d_eig2 <- multiomicCCA::compute_laplacian(protein_frnn$d_g, k_max = k_max, rowname_vec = cell_name, 
                                          colname_vec = paste0("dlap_", 1:k_max), verbose = F)
e_eig2 <- multiomicCCA::compute_laplacian(protein_frnn$e_g, k_max = k_max, rowname_vec = cell_name, 
                                          colname_vec = paste0("elap_", 1:k_max), verbose = F)

set.seed(10)
protein_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = F, data_2 = T,
                                                     c_g = protein_frnn$c_g, d_g = protein_frnn$d_g, 
                                                     only_embedding = T, verbose = T, 
                                                     sampling_type = "adaptive_gaussian",
                                                     keep_nn = T)

save(date_of_run, session_info, protein_frnn, c_eig2, d_eig2, e_eig2, protein_embeddings,
     file = "../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca3.RData")

rm(list = "protein_embeddings"); gc(T)

##########################

set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, g_1 = rna_frnn$c_g,
                                         g_2 = protein_frnn$c_g, nn = nn, 
                                         verbose = 2)
combined_g <- SeuratObject::as.Graph(combined_g)
set.seed(10)
combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings

save(date_of_run, session_info, combined_common_umap,
     file = "../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca4.RData")

rm(list = c("rna_frnn", "protein_frnn", "combined_g")); gc(T)

set.seed(10)
both_embeddings <- multiomicCCA::plot_embeddings(dcca_res, membership_vec = membership_vec,
                                                 data_1 = T, data_2 = T, add_noise = F,
                                                 only_embedding = T)

save(date_of_run, session_info, both_embeddings,
     file = "../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca5.RData")

rm(list = c("dcca_res", "both_embeddings")); gc(T)
#####################################

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, rank_c = min(rank_1, rank_2), 
                                                verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
mat_2_denoised <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2

rm(list = c("dcca_res", "dcca_decomp")); gc(T)

p1 <- ncol(mat_1_denoised)
gene_smoothed <- lapply(1:p1, function(j){
  print(j)
  
  c_res <- compute_smooth_signal(mat_1_denoised[,j], c_eig)
  d_res <- compute_smooth_signal(mat_1_denoised[,j], d_eig)
  e_res <- compute_smooth_signal(mat_1_denoised[,j], e_eig)
  
  list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
       d_variance = d_res$variance, d_r2 = d_res$r_squared,
       e_variance = e_res$variance, e_r2 = e_res$r_squared)
})

rm(list = c("mat_1_denoised", "c_eig", "d_eig", "e_eig")); gc(T)

p2 <- ncol(mat_2_denoised)
protein_smoothed <- lapply(1:p2, function(j){
  print(j)
  
  c_res <- compute_smooth_signal(mat_2_denoised[,j], c_eig2)
  d_res <- compute_smooth_signal(mat_2_denoised[,j], d_eig2)
  e_res <- compute_smooth_signal(mat_2_denoised[,j], e_eig2)
  
  list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
       d_variance = d_res$variance, d_r2 = d_res$r_squared,
       e_variance = e_res$variance, e_r2 = e_res$r_squared)
})

save(date_of_run, session_info, gene_smoothed, protein_smoothed,
     file = "../../../../out/Writeup14b/Writeup14b_citeseq_pbmc228_dcca6.RData")


