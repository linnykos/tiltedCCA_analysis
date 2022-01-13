rm(list=ls())
set.seed(10)

mat1 <- matrix(rnorm(100), 10, 10)
mat2 <- matrix(rnorm(100), 10, 10)
u1 <- svd(mat1)$u[,1:2]
u2 <- svd(mat2)$u[,1:5]

svd_tmp <- svd(t(u1) %*% u2)
svd_tmp$d
