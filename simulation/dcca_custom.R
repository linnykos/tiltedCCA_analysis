dcca_custom <- function(dcca_res){
  
  common_1 <- multiomicCCA:::.extract_matrix_helper(common_score = dcca_res$common_score, 
                                                    distinct_score = dcca_res$distinct_score_1,
                                                    svd_e = dcca_res$svd_1, 
                                                    common_bool = T, 
                                                    distinct_bool = F,
                                                    center = F, 
                                                    renormalize = F)
  distinct_1 <- multiomicCCA:::.extract_matrix_helper(common_score = dcca_res$common_score, 
                                                      distinct_score = dcca_res$distinct_score_1,
                                                      svd_e = dcca_res$svd_1, 
                                                      common_bool = F, 
                                                      distinct_bool = T,
                                                      center = F, 
                                                      renormalize = F)
  common_2 <- multiomicCCA:::.extract_matrix_helper(common_score = dcca_res$common_score, 
                                                    distinct_score = dcca_res$distinct_score_2,
                                                    svd_e = dcca_res$svd_2, 
                                                    common_bool = T, 
                                                    distinct_bool = F,
                                                    center = F, 
                                                    renormalize = F)
  distinct_2 <- multiomicCCA:::.extract_matrix_helper(common_score = dcca_res$common_score, 
                                                      distinct_score = dcca_res$distinct_score_2,
                                                      svd_e = dcca_res$svd_2, 
                                                      common_bool = F, 
                                                      distinct_bool = T,
                                                      center = F, 
                                                      renormalize = F)
  
  ##########
  common <- .pca_custom(cbind(common_1, common_2))
  distinct_1 <- .pca_custom(distinct_1)
  distinct_2 <- .pca_custom(distinct_2)
  
  list(common = common,
       distinct_1 = distinct_1,
       distinct_2 = distinct_2)
}

.pca_custom <- function(mat){
  mat <- scale(mat, center = T, scale = F)
  svd_res <- svd(mat)
  multiomicCCA:::.mult_mat_vec(svd_res$u[,1:2], svd_res$d[1:2])
}
