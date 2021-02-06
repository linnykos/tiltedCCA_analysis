rm(list=ls())
vec1 <- c(1,5); vec2 <- c(4,2)
vec1 <- vec1/.l2norm(vec1)*20
vec2 <- vec2/.l2norm(vec2)*10

.decomposition_2d(vec1, vec2, plotting = T)

##########################

rm(list=ls())
vec1 <- c(1,5); vec2 <- c(10,1)
vec1 <- vec1/.l2norm(vec1)*3
vec2 <- vec2/.l2norm(vec2)

.decomposition_2d(vec1, vec2, plotting = T)

##########################

rm(list=ls())
vec1 <- c(1,10); vec2 <- c(10,1)
vec1 <- vec1/.l2norm(vec1)
vec2 <- vec2/.l2norm(vec2)

.decomposition_2d(vec1, vec2, plotting = T)

##########################

rm(list=ls())
vec1 <- c(1,1); vec2 <- c(1.1,0.9)
vec1 <- vec1/.l2norm(vec1)*3
vec2 <- vec2/.l2norm(vec2)

.decomposition_2d(vec1, vec2, plotting = T)
