rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_common <- rna_common[cell_idx[order_vec],]
atac_pred <- atac_pred[cell_idx[order_vec],]

p <- ncol(rna_common); n <- nrow(rna_common)
rna_predicted_mat <- sapply(1:p, function(j){
  if(j > 10 & j %% floor(p/10) == 0) cat('*')
  
  tmp_df <- data.frame(y = rna_common[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
atac_predicted_mat <- sapply(1:p, function(j){
  if(j > 10 & j %% floor(p/10) == 0) cat('*')
  
  tmp_df <- data.frame(y = atac_pred[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})

save(greenleaf, multiSVD_obj, 
     rna_common, atac_pred, 
     rna_predicted_mat, atac_predicted_mat,
     date_of_run, session_info,
     file = "../../../out/main/10x_greenleaf_timeseries.RData")
