rm(list=ls())
x=1
set.seed(x)
vec1 <- abs(rnorm(2)); vec1 <- vec1/.l2norm(vec1)
vec2 <- abs(rnorm(2)); vec2 <- vec2/.l2norm(vec2)
circle <- .construct_circle(vec1, vec2)
tmp <- .rightmost_vector(vec1, vec2)
vec1 <- tmp$vec_left; vec2 <- tmp$vec_right
lower_radian <- .find_radian(circle, tmp$vec_right)
upper_radian <- .find_radian(circle, tmp$vec_left)
stopifnot(lower_radian < upper_radian)
distinct_perc_2 <- runif(1)
distinct_perc_2

res <- .binary_search_radian(circle, lower_radian, upper_radian, 
                             distinct_perc_2, max_iter = 50, tol = 1e-6)
common_vec <- .position_from_circle(circle, res)

distinct1 <- vec1 - common_vec
distinct2 <- vec2 - common_vec
ratio <- .l2norm(distinct2)/(.l2norm(distinct1) + .l2norm(distinct2))
ratio

plot(NA, xlim = c(0,1), ylim = c(0,1), asp = T)
lines(c(0,vec1[1]), c(0,vec1[2]))
lines(c(0,vec2[1]), c(0,vec2[2]))
lines(c(0,common_vec[1]), c(0, common_vec[2]), col = "green")

##########################

lower <- lower_radian; upper <- upper_radian
right_vec <- .position_from_circle(circle, lower_radian)
left_vec <- .position_from_circle(circle, upper_radian)
iter <- 1; prev_mid <- NA
tol <- 1e-6; max_iter = 10

while(iter < max_iter){
  mid <- (lower+upper)/2
  print(paste0("mid: ", round(mid,2)))
  common_vec <- .position_from_circle(circle, mid)
  left_distinct <- left_vec - common_vec
  right_distinct <- right_vec - common_vec
  ratio <- .l2norm(left_distinct)/(.l2norm(left_distinct)+.l2norm(right_distinct))
  print(paste0("ratio: ", round(ratio,2), ", target: ", round(distinct_perc_2, 2)))
  if(abs(ratio - distinct_perc_2) <= tol) break()
  
  if(ratio > distinct_perc_2){
    # need to make left_distinct smaller, so move radian right (i.e., smaller radian)
    lower <- mid
  } else {
    upper <- mid
  }
  if(!is.na(prev_mid) && abs(mid - prev_mid) < tol) break()
  
  prev_mid <- mid
  iter <- iter+1
}
