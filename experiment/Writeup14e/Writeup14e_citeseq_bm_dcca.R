rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)

#########################

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

##################

rank_1 <- 30; rank_2 <- 18; nn <- 30
svd_1 <- multiomicCCA:::.svd_truncated(mat_1b, K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = T, sd_vec = T,
                                       symmetric = F)
tmp_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
snn_mat_1 <- multiomicCCA:::.form_snn_mat(bool_intersect = T,
                                          mat = tmp_1, 
                                          num_neigh = nn)
metacell_clustering_1 <- lapply(1:nrow(snn_mat_1), function(i){
  which(snn_mat_1[i,] != 0)
})

svd_2 <- multiomicCCA:::.svd_truncated(mat_2b, K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
tmp_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u[,-1], svd_2$d[-1])
snn_mat_2 <- multiomicCCA:::.form_snn_mat(bool_intersect = T,
                                          mat = tmp_2, 
                                          num_neigh = nn)
metacell_clustering_2 <- lapply(1:nrow(snn_mat_2), function(i){
  which(snn_mat_2[i,] != 0)
})

###########################

set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, 
                                      dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      center_1 = T, center_2 = T,
                                      scale_1 = T, scale_2 = T,
                                      num_neigh = nn, 
                                      metacell_clustering_1 = metacell_clustering_1,
                                      metacell_clustering_2 = metacell_clustering_2,
                                      fix_tilt_perc = F, verbose = T)
dcca_res$cca_obj
dcca_res$df_percentage
dcca_res$tilt_perc

dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, verbose = T)
dcca_res2$tilt_perc

save(bm, dcca_res, dcca_res2, metacell_clustering_1, metacell_clustering_2,
     rank_1, rank_2, nn, date_of_run, session_info,
     file = "../../../../out/Writeup14e/Writeup14e_citeseq_bm_dcca.RData")
