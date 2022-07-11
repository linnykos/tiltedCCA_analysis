rm(list=ls())
set.seed(10)

n_each <- 100
true_membership_vec <- rep(1:3, each = n_each)
clustering_1 <- factor(c(rep(1, 2*n_each), rep(2, n_each)))
clustering_2 <- factor(c(rep(1, n_each), rep(2, n_each), rep(1, n_each)))
mat_1 <- do.call(rbind, lapply(1:3, function(i){
  if(i %in% c(1,2)){
    MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(3,0), Sigma = diag(2)) 
  }
}))

mat_2 <- do.call(rbind, lapply(1:3, function(i){
  if(i %in% c(1,3)){
    MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(2.5,0), Sigma = diag(2)) 
  }
}))

mat_1 <- scale(mat_1, center = T, scale = F)
mat_2 <- scale(mat_2, center = T, scale = F)
svd_1 <- svd(mat_1)
svd_2 <- svd(mat_2)

p_1 <- 40; p_2 <- 40
svd_v_1 <- generate_random_orthogonal(p_1, 2)
svd_v_2 <- generate_random_orthogonal(p_2, 2)

mat_1 <- tcrossprod(.mult_mat_vec(svd_1$u, svd_1$d), svd_v_1)
mat_2 <- tcrossprod(.mult_mat_vec(svd_2$u, svd_2$d), svd_v_2)
mat_1 <- mat_1 + rnorm(prod(dim(mat_1)), sd = 0.5)
mat_2 <- mat_2 + rnorm(prod(dim(mat_2)), sd = 0.5)

rownames(mat_1) <- paste0("n", 1:nrow(mat_1))
rownames(mat_2) <- paste0("n", 1:nrow(mat_2))
colnames(mat_1) <- paste0("g", 1:ncol(mat_1))
colnames(mat_2) <- paste0("p", 1:ncol(mat_2))

##############

n <- nrow(mat_1)
large_clustering_1 <- clustering_1
large_clustering_2 <- clustering_2
multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                dims_1 = 1:3, dims_2 = 1:3,
                                center_1 = F, center_2 = F,
                                normalize_row = T,
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
                             num_neigh = 10,
                             bool_cosine = T,
                             bool_intersect = T,
                             min_deg = 1)

multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj)
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1)
  
###############

set.seed(10)
col_vec <- colorRampPalette(c("white", 
                              rgb(98, 147, 194, maxColorValue = 255),
                              rgb(65, 88, 163, maxColorValue = 255)))(10)
tmp <- mat_1[1:13,1:10]
png("../../out/figures/main/toy_simulation_rna-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = 12/9, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

col_vec <- colorRampPalette(c("white", 
                              rgb(171, 213, 113, maxColorValue = 255),
                              rgb(82, 149, 59, maxColorValue = 255)))(10)
tmp <- mat_2[1:7,1:10]
png("../../out/figures/main/toy_simulation_adt-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = 6/9, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

col_vec <- colorRampPalette(c("white", 
                              rgb(219, 154, 237, maxColorValue = 255),
                              rgb(191, 74, 223, maxColorValue = 255)))(10)
tmp <- mat_2[1:13,1:20]
png("../../out/figures/main/toy_simulation_atac-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = 12/19, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

set.seed(10)
col_vec <- colorRampPalette(c("white", 
                              rgb(246, 224, 157, maxColorValue = 255),
                              rgb(236, 170, 16, maxColorValue = 255)))(10)
tmp <- multiSVD_obj$common_mat_1[1:13,1:10]
png("../../out/figures/main/toy_simulation_rna-commonMat-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = 12/9, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

set.seed(10)
col_vec <- colorRampPalette(c("white", 
                              rgb(98, 147, 194, maxColorValue = 255),
                              rgb(65, 88, 163, maxColorValue = 255)))(10)
tmp <- multiSVD_obj$common_mat_1[1:13,1:10]
png("../../out/figures/main/toy_simulation_rna-commonMat-heatmap2.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = 12/9, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

set.seed(10)
col_vec <- colorRampPalette(c("white", 
                              rgb(98, 147, 194, maxColorValue = 255),
                              rgb(65, 88, 163, maxColorValue = 255)))(10)
tmp <- multiSVD_obj$distinct_mat_1[1:13,1:10]
png("../../out/figures/main/toy_simulation_rna-distinctMat-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = 12/9, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

set.seed(10)
col_vec <- colorRampPalette(c("white", 
                              rgb(246, 224, 157, maxColorValue = 255),
                              rgb(236, 170, 16, maxColorValue = 255)))(10)
tmp <- multiSVD_obj$common_mat_2[1:7,1:10]
png("../../out/figures/main/toy_simulation_adt-commonMat-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = 6/9, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

set.seed(10)
col_vec <- colorRampPalette(c("white", 
                              rgb(171, 213, 113, maxColorValue = 255),
                              rgb(82, 149, 59, maxColorValue = 255)))(10)
tmp <- multiSVD_obj$distinct_mat_2[1:7,1:10]
png("../../out/figures/main/toy_simulation_adt-distinctMat-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(tmp), asp = 6/9, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

col_vec <- colorRampPalette(c("white", 
                              rgb(246, 224, 157, maxColorValue = 255),
                              rgb(236, 170, 16, maxColorValue = 255)))(10)
png("../../out/figures/main/toy_simulation_commonScore-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(multiSVD_obj$tcca_obj$common_score[1:10,], asp = 2/9, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()


col_vec <- colorRampPalette(c("white", 
                              rgb(98, 147, 194, maxColorValue = 255),
                              rgb(65, 88, 163, maxColorValue = 255)))(10)
png("../../out/figures/main/toy_simulation_rna-distinctScore-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(multiSVD_obj$tcca_obj$distinct_score_1[1:10,], asp = 2/9, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

col_vec <- colorRampPalette(c("white", 
                              rgb(171, 213, 113, maxColorValue = 255),
                              rgb(82, 149, 59, maxColorValue = 255)))(10)
png("../../out/figures/main/toy_simulation_adt-distinctScore-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(multiSVD_obj$tcca_obj$distinct_score_2[1:10,], asp = 2/9, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

col_vec <- colorRampPalette(c("white", 
                              rgb(0.65, 0.65, 0.65)))(10)
png("../../out/figures/main/toy_simulation_rna-coefficient-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(mat_1[1:13,11:13]), asp = 9/2, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

col_vec <- colorRampPalette(c("white", 
                              rgb(0.65, 0.65, 0.65)))(10)
png("../../out/figures/main/toy_simulation_adt-coefficient-heatmap.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
image(t(mat_2[1:7,11:13]), asp = 6/2, col = col_vec, 
      xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

#############

png("../../out/figures/main/toy_simulation_emptyAxis.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(1, 1, 1, 1))
plot(NA, xlim = c(0,1), ylim = c(0,1), bty = "n", xaxt = "n", yaxt = "n",
     xlab = "", ylab = "", asp = T)
axis(1, at = seq(0,1,length.out=5), labels = rep("",5), 
     lwd = 3, lwd.ticks = 3)
axis(2, at = seq(0,1,length.out=5), labels = rep("",5), 
     lwd = 3, lwd.ticks = 3)
graphics.off()

plot_decomposition_2d(vec1 = multiSVD_obj$cca_obj$score_1[,1],
                      vec2 = multiSVD_obj$cca_obj$score_2[,1],
                      common_vec = multiSVD_obj$tcca_obj$common_score[,1],
                      xlim = c(0,2),
                      ylim = c(0,2))

image(as.matrix(multiSVD_obj$snn_list$snn_1), asp = T)
image(as.matrix(multiSVD_obj$snn_list$snn_2), asp = T)

#######################

set.seed(10)

n_each <- 5
clustering_1 <- factor(c(rep(1, 2*n_each), rep(2, n_each)))
clustering_2 <- factor(c(rep(1, 3), rep(2, 5), rep(1, 2)))
mat_1 <- do.call(rbind, lapply(1:2, function(i){
  if(i == 1){
    MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(2,0), Sigma = diag(2)) 
  }
}))

mat_2 <- do.call(rbind, lapply(1:2, function(i){
  if(i == 1){
    MASS::mvrnorm(n = n_each, mu = c(0,0), Sigma = diag(2)) 
  } else {
    MASS::mvrnorm(n = n_each, mu = c(0.5,0), Sigma = diag(2)) 
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
multiSVD_obj2 <- compute_snns(input_obj = multiSVD_obj,
                             latent_k = 2,
                             num_neigh = 2,
                             bool_cosine = F,
                             bool_intersect = T,
                             min_deg = 1)

spacing <- 1/((10-1)*2)
xseq <- seq(-spacing, 1+spacing, by = 2*spacing)

png("../../out/figures/main/toy_simulation_rna-SNN.png", 
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

png("../../out/figures/main/toy_simulation_adt-SNN.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(as.matrix(multiSVD_obj2$snn_list$snn_2)), asp = T,
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
tmp2 <- as.matrix(multiSVD_obj2$snn_list$snn_2)
tmp3 <- tmp1 + tmp2; tmp3[tmp3 > 0] <- 1
png("../../out/figures/main/toy_simulation_Common-SNN.png", 
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

tmp4 <- tmp3
set.seed(10)
samp_idx <- cbind(sample(1:n), sample(1:n))
zz <- samp_idx[,1] != samp_idx[,2]
samp_idx <- samp_idx[zz,]
for(k in 1:nrow(samp_idx)){
  i <- samp_idx[k,1]; j <- samp_idx[k,2]
  tmp4[i,j] <- -tmp4[i,j]+1
  tmp4[j,i] <- -tmp4[j,i]+1
}
png("../../out/figures/main/toy_simulation_Attmpted-SNN.png", 
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
