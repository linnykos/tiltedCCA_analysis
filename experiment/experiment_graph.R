mat = matrix(c(0,1,0,0,0, 
               1,0,0,1,0,
               0,0,0,1,1,
               0,1,1,0,1,
               0,0,1,1,0), 5, 5)
eigen(mat)

deg <- rowSums(mat)
lap <- diag(deg) - mat
eigen(lap)
lap

lap <- diag(deg^(-.5)) %*% mat %*% diag(deg^(-.5))
round(lap,2)
zz <- eigen(lap)
zz
yy <- diag(deg^(-1/2)) %*% zz$vectors
for(j in 1:ncol(yy)){
  yy[,j] <- yy[,j]/sqrt(sum(yy[,j]^2))
}
yy


lap <- diag(5) - diag(deg^(-.5)) %*% mat %*% diag(deg^(-.5))
round(lap,2)
eigen(lap)
eigen(lap)$vectors[,5]/sqrt(deg)

lap <- diag(deg^(-1))%*%mat
round(lap,2)
zz <- eigen(lap, symmetric = F)
zz
lap %*% zz$vectors[,3]
zz$values[3] * zz$vectors[,3]


lap <- diag(5) - diag(deg^(-1))%*%mat
round(lap,2)
eigen(lap, symmetric = F)


lap <-  diag(deg^(-1))%*%mat %*% diag(deg^(-1))
deg2 <- rowSums(lap)
lap <- diag(deg2^(-1))%*%lap
zz <- eigen(lap)
