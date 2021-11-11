rm(list=ls())

load("../../../../out/Writeup14f/Writeup14f_SNU_preprocessed.RData")
date_of_run <- Sys.time(); session_info <- devtools::session_info()

Seurat::DefaultAssay(SNU) <- "atac"
mat_1 <- Matrix::t(SNU[["atac"]]@data[Seurat::VariableFeatures(object = SNU),])
Seurat::DefaultAssay(SNU) <- "cna"
mat_2 <- Matrix::t(SNU[["cna"]]@data)

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

rank_1 <- 50; rank_2 <- 40; nn <- 30
set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1b, mat_2b, 
                                      dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      num_neigh = nn, 
                                      metacell_clustering_1 = NA,
                                      metacell_clustering_2 = NA,
                                      fix_tilt_perc = F, verbose = T)
dcca_res$cca_obj
dcca_res$df_percentage
dcca_res$tilt_perc

dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, verbose = T)
dcca_res2$tilt_perc

save(SNU, dcca_res, dcca_res2, 
     rank_1, rank_2, nn, date_of_run, session_info,
     file = "../../../../out/Writeup14f/Writeup14f_SNU601_dcca.RData")
