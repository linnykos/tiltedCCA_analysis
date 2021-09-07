rm(list=ls())
n <- 5; K <- 2
U <- matrix(runif(n*K), nrow = n, ncol = K)
V <- matrix(runif(n*K), nrow = n, ncol = K)
for(i in 1:5){
  U[i,] <- U[i,]/sum(U[i,])
}
for(i in 1:2){
  V[,i] <- V[,i]/sum(V[,i])
}
U <- round(U,2); V <- round(V,2)

zz1 <- U%*%t(V); rowSums(zz1)
zz2 <- t(V)%*%U; rowSums(zz2)
