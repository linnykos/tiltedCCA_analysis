m <- 2
P = matrix(c(1-1/(2*m), 1/(2*m), 0, 0, 0,
             1/2, 0, 1/2, 0, 0,
             0, 1/2, 0, 1/2, 0,
             0, 0, 1/2, 0, 1/2,
             0, 0, 0, 1/(2*m), 1-1/(2*m)), 5, 5, byrow = T)
res <- eigen(P)
Z <- 2*m+3
stationary <- c(m,1,1,1,m)/Z
.l2norm(stationary)
stationary/.l2norm(stationary)

res

plot(res$vectors[,1], ylim = range(res$vectors))
res$vectors %*% diag(res$values) %*% solve(res$vectors)
round(crossprod(res$vectors),2) #notice: not orthogonal
# the stationary distribution is hidden in solve(res$vectors)[1,]

iter <- 100
P2 <- P
for(i in 1:iter){
  P2 <- P2%*%P
}
plot(P2[1,], ylim = range(res$vectors))

.diag_mat <- function(x){
  if(length(x) == 1) return(matrix(x,1,1))
  diag(x)
}

for(i in 1:5){
  print(i)
  left_vec <- res$vectors[,1:i,drop = F]
  val <- res$values[1:i]
  right_vec <- MASS::ginv(left_vec)
  approx_mat <- left_vec %*% .diag_mat(val) %*% right_vec
  print(.l2norm(P - approx_mat))
}
