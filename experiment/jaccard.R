trials <- 10000
set.seed(10)
d1 <- 1; d2 <- 20; c <- 5
mat <- matrix(runif(2*trials), ncol = 2, nrow = trials)
mat <- rbind(mat, c(1,1), c(0,0))
vec <- apply(mat, 1, function(x){
  alpha <- x[1]; beta <- x[2]
  min(c((alpha*d1+c)/(d1+beta*d2+c), (beta*d2+c)/(alpha*d1+d2+c)))
})
idx <- which.max(vec)
max(vec)
mat[idx,]

###########3

set.seed(10)
d1 <- 10; d2 <- 20
mat <- matrix(runif(2*trials), ncol = 2, nrow = trials)
mat <- rbind(mat, c(1,1), c(1, 0.65))
vec <- apply(mat, 1, function(x){
  alpha <- x[1]; beta <- x[2]
  min(c(alpha*d1/(d1+beta*d2), beta*d2/(alpha*d1+d2)))
})
idx <- which.max(vec)
max(vec)
mat[idx,]