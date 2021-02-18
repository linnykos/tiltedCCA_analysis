mat_1 <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
mat_2 <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2
n <- nrow(mat_1)

coef_mat_1 <- crossprod(dcca_res$score_1, mat_1)/n
coef_mat_2 <- crossprod(dcca_res$score_2, mat_2)/n

zz <- dcca_res$distinct_score_1[,1,drop = F] %*% coef_mat_1[1,,drop = F]
zz <- multiomicCCA:::.svd_truncated(zz, K = 1)
zz$d
multiomicCCA:::.l2norm(dcca_res$distinct_score_1[,1])*multiomicCCA:::.l2norm(coef_mat_1[1,])

zz <- dcca_res$common_score[,1,drop = F] %*% coef_mat_1[1,,drop = F]
zz <- multiomicCCA:::.svd_truncated(zz, K = 1)
zz$d
multiomicCCA:::.l2norm(dcca_res$common_score[,1])*multiomicCCA:::.l2norm(coef_mat_1[1,])
