rm(list=ls())
vec1 <- c(1,5); vec2 <- c(4,2)
vec1 <- vec1/.l2norm(vec1)*20
vec2 <- vec2/.l2norm(vec2)*10

png("../../out/simulation/Writeup11/Writeup11_dmca_decomposition.png", height = 900, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3))
.decomposition_2d(vec1, vec2, plotting = T)

##########################

vec1 <- c(1,5); vec2 <- c(10,1)
vec1 <- vec1/.l2norm(vec1)*3
vec2 <- vec2/.l2norm(vec2)

.decomposition_2d(vec1, vec2, plotting = T)

##########################

# rm(list=ls())
# vec1 <- c(1,10); vec2 <- c(10,1)
# vec1 <- vec1/.l2norm(vec1)
# vec2 <- vec2/.l2norm(vec2)
# 
# .decomposition_2d(vec1, vec2, plotting = T)

##########################

vec1 <- c(1,1); vec2 <- c(1.1,0.9)
vec1 <- vec1/.l2norm(vec1)*3
vec2 <- vec2/.l2norm(vec2)

.decomposition_2d(vec1, vec2, plotting = T)
graphics.off()
