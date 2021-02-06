rm(list=ls())
x <- 2
set.seed(x)
low_rank1 <- matrix(rnorm(50), nrow = 25, ncol = 2)
mat1 <- tcrossprod(low_rank1)
svd_1 <- .svd_truncated(mat1, K = 2)
mat1 <- .recover_mat_from_svd(svd_1)

low_rank2 <- matrix(rnorm(50), nrow = 25, ncol = 2)
mat2 <- tcrossprod(low_rank2)
svd_2 <- .svd_truncated(mat2, K = 2)
mat2 <- .recover_mat_from_svd(svd_2)

mca_res <- .mca(svd_1, svd_2, rank = 2)
# res <- .mca_common_score(svd_1, svd_2, mca_res, verbose = F)
###################

n <- nrow(svd_1$u); K <- ncol(mca_res$u)
score_1 <- .mult_mat_vec(svd_1$u, svd_1$d) %*% crossprod(svd_1$v, mca_res$u)
score_2 <- .mult_mat_vec(svd_2$u, svd_2$d) %*% crossprod(svd_2$v, mca_res$v)

common_mat_1 <- matrix(NA, nrow = n, ncol = K)
common_mat_2 <- matrix(NA, nrow = n, ncol = K)
j <- 1
basis_res <- .representation_2d(score_1[,j], score_2[,j])
tmp_decomp <- .decomposition_2d(basis_res$rep1, basis_res$rep2, plotting = F)
