pca_combine <- function(svd_1,
                        svd_2,
                        dims_1,
                        dims_2){
  dimred1 <- multiomicCCA:::.mult_mat_vec(svd_1$u[,dims_1,drop = F], svd_1$d[dims_1]/svd_1$d[1])
  dimred2 <- multiomicCCA:::.mult_mat_vec(svd_2$u[,dims_2,drop = F], svd_2$d[dims_2]/svd_1$d[2])
  
  dimred_all <- cbind(dimred1, dimred2)
  dimred_all <- scale(dimred_all, scale = F, center = T)
  svd_res <- svd(dimred_all)
  tmp <- multiomicCCA:::.mult_mat_vec(svd_res$u[,1:2,drop = F], svd_res$d[1:2])
  tmp
}