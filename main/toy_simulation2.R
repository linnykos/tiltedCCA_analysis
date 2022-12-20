rm(list=ls())
set.seed(10)

n_each <- 4
clustering_1 <- factor(c(rep(1, 2*n_each), rep(2, n_each)))
clustering_2 <- factor(c(rep(1, 3), rep(2, 5), rep(1, 2)))
mat_1 <- do.call(rbind, lapply(1:3, function(i){
  if(i == 1){
    MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
  } else if (i == 2){
    MASS::mvrnorm(n = n_each, mu = c(0,5), Sigma = diag(2)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(5,0), Sigma = diag(2)) 
  }
}))

mat_2 <- do.call(rbind, lapply(1:3, function(i){
  if(i == 1){
    MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(5,0), Sigma = diag(2)) 
  }
}))

n <- nrow(mat_1)
large_clustering_1 <- clustering_1
large_clustering_2 <- clustering_2
multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                dims_1 = 1:2, dims_2 = 1:2,
                                center_1 = F, center_2 = F,
                                normalize_row = F,
                                normalize_singular_value = F,
                                recenter_1 = F, recenter_2 = F,
                                rescale_1 = F, rescale_2 = F,
                                scale_1 = F, scale_2 = F)
multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                               large_clustering_1 = NULL, 
                               large_clustering_2 = NULL,
                               num_metacells = NULL)
multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                             latent_k = 2,
                             num_neigh = 3,
                             bool_cosine = F,
                             bool_intersect = T,
                             min_deg = 1)

##################################


spacing <- 1/((nrow(mat_1)-1)*2)
xseq <- seq(-spacing, 1+spacing, by = 2*spacing)

png("../../out/main/toy_SNN_1.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(as.matrix(multiSVD_obj$snn_list$snn_1)), asp = T,
      col = c(rgb(0.9, 0.9, 0.9), 2),
      breaks = c(-.5, .5, 1.5),
      xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:length(xseq)){
  lines(rep(xseq[i], 2), range(xseq), lwd = 2)
}
for(i in 1:length(xseq)){
  lines(range(xseq), rep(xseq[i], 2), lwd = 2)
}
graphics.off()

png("../../out/main/toy_SNN_2.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(as.matrix(multiSVD_obj$snn_list$snn_2)), asp = T,
      col = c(rgb(0.9, 0.9, 0.9), 2),
      breaks = c(-.5, .5, 1.5),
      xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:length(xseq)){
  lines(rep(xseq[i], 2), range(xseq), lwd = 2)
}
for(i in 1:length(xseq)){
  lines(range(xseq), rep(xseq[i], 2), lwd = 2)
}
graphics.off()

tmp1 <- as.matrix(multiSVD_obj$snn_list$snn_1)
tmp2 <- as.matrix(multiSVD_obj$snn_list$snn_2)
tmp3 <- tmp1 + tmp2; tmp3[tmp3 > 0] <- 1
png("../../out/main/toy_common_SNN.png",  
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(tmp3), asp = T,
      col = c(rgb(0.9, 0.9, 0.9), 2),
      breaks = c(-.5, .5, 1.5),
      xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:length(xseq)){
  lines(rep(xseq[i], 2), range(xseq), lwd = 2)
}
for(i in 1:length(xseq)){
  lines(range(xseq), rep(xseq[i], 2), lwd = 2)
}
graphics.off()

tmp4 <- as.matrix(multiSVD_obj$snn_list$common_snn)
set.seed(10)
samp_idx <- cbind(sample(1:n), sample(1:n))
zz <- samp_idx[,1] != samp_idx[,2]
samp_idx <- samp_idx[zz,]
for(k in 1:nrow(samp_idx)){
  i <- samp_idx[k,1]; j <- samp_idx[k,2]
  tmp4[i,j] <- -tmp4[i,j]+1
  tmp4[j,i] <- -tmp4[j,i]+1
}
for(k in 1:nrow(tmp4)){
  idx <- which(tmp4[i,] != 0)
  if(length(idx) > 4){
    tmp4[i,] <- 0
    tmp4[i,sample(idx,3)] <- 1
  }
}
tmp4 <- tmp4 + t(tmp4)
tmp4[tmp4 > 0.5] <- 1
png("../../out/main/toy_attempted_SNN.png",   
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(tmp4), asp = T,
      col = c(rgb(0.9, 0.9, 0.9), 2),
      breaks = c(-.5, .5, 1.5),
      xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:length(xseq)){
  lines(rep(xseq[i], 2), range(xseq), lwd = 2)
}
for(i in 1:length(xseq)){
  lines(range(xseq), rep(xseq[i], 2), lwd = 2)
}
graphics.off()