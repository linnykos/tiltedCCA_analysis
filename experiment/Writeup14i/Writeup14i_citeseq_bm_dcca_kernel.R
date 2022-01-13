rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
library(multiomicCCA)
source("kernel_matrix.R")
source("metacell.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(bm)

Seurat::DefaultAssay(bm) <- "RNA"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30)
bm <- Seurat::FindClusters(bm, resolution = 0.25)

Seurat::DefaultAssay(bm) <- "ADT"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:18, reduction = "apca")
bm <- Seurat::FindClusters(bm, resolution = 0.25)

###############

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)

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

########################

clustering_1 <- lapply(levels(bm$RNA_snn_res.0.25), function(x){
  which(bm$RNA_snn_res.0.25 == x)
})
clustering_1 <- clustering_1[sapply(clustering_1, length) != 0]
clustering_2 <- lapply(levels(bm$ADT_snn_res.0.25), function(x){
  which(bm$ADT_snn_res.0.25 == x)
})
clustering_2 <- clustering_2[sapply(clustering_2, length) != 0]

metacell_clustering_1 <- form_subclusters(mat = mat_1b, 
                                          clustering = clustering_1, 
                                          target_k = 80)
metacell_clustering_2 <- form_subclusters(mat = mat_2b, 
                                          clustering = clustering_2, 
                                          target_k = 80)
tmp <- intersect_metacells(metacell_clustering_1 = metacell_clustering_1,
                           metacell_clustering_2 = metacell_clustering_2,
                           n = ncol(bm))
metacell_clustering <- tmp$metacell_clustering
clustering_hierarchy_1 <- tmp$clustering_hierarchy_1
clustering_hierarchy_2 <- tmp$clustering_hierarchy_2
print(paste0(length(metacell_clustering), " number of metacells"))

#########################

dist_mat_1 <- form_dist_matrix(mat = mat_1, K = 30, 
                               metacell_clustering = metacell_clustering)

dist_mat_2 <- form_dist_matrix(mat = mat_2, K = 18, 
                               metacell_clustering = metacell_clustering)

tmp <- compute_min_embedding(dist_mat_1 = dist_mat_1, 
                                       dist_mat_2 = dist_mat_2)
min_dist_mat <- tmp$dist_mat
target_embedding <- tmp$dimred

#########################

rank_1 <- 30; rank_2 <- 18; nn <- 30
set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1b, mat_2b, 
                                      dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      discretization_gridsize = 11,
                                      num_neigh = nn, 
                                      metacell_clustering = metacell_clustering,
                                      fix_tilt_perc = F, 
                                      enforce_boundary = F,
                                      target_embedding = target_embedding,
                                      verbose = T)
dcca_res$cca_obj
dcca_res$df_percentage
dcca_res$tilt_perc

dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, verbose = T)
dcca_res2$tilt_perc

save(bm, dcca_res, dcca_res2, 
     rank_1, rank_2, nn, date_of_run, session_info,
     dist_mat_1, dist_mat_2,
     metacell_clustering_1, 
     metacell_clustering_2,
     metacell_clustering,
     min_dist_mat, target_embedding,
     clustering_hierarchy_1,
     clustering_hierarchy_2,
     file = "../../../../out/Writeup14i/Writeup14i_citeseq_bm_dcca.RData")


