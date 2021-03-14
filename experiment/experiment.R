rm(list=ls())
trials <- 100
rank_vec <- sapply(1:trials, function(x){
  set.seed(x)
  B <- matrix(rnorm(9), 3, 3)
  B <- scale(B, center = T, scale = F)
  Matrix::rankMatrix(B)
})