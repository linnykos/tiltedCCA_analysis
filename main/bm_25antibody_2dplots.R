rm(list=ls())
load("../../../out/main/citeseq_bm25.RData")
source("bm_25antibody_colorPalette.R")

celltype_vec <- bm$celltype.l2
col_vec <- sapply(celltype_vec, function(celltype){
  col_palette[which(names(col_palette) == celltype)]
})
col_vec2 <- scales::col2hcl(col_vec, l = 100)

x_vec <- multiSVD_obj$cca_obj$score_1[,1]
y_vec <- multiSVD_obj$cca_obj$score_2[,1]
diag_vec <- rep(1/sqrt(2),2)
projection_vec <- cbind(x_vec, y_vec) %*% diag_vec %*% t(diag_vec)
png("../../../out/figures/main/citeseq_bm25_2d_cca.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(1, 1, 1, 1))
plot(x_vec, y_vec,
     asp = T, xaxt = "n", yaxt = "n", bty = "n",
     col = col_vec2, pch = 16)
points(projection_vec[,1], projection_vec[,2], col = "white", cex = 2, pch = 16)
points(projection_vec[,1], projection_vec[,2], col = col_vec, pch = 16)
axis(side = 1, labels = F)
axis(side = 2, labels = F)
graphics.off()

png("../../../out/figures/main/citeseq_bm25_2d_cca_cleaned.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(1, 1, 1, 1))
plot(x_vec, y_vec,
     asp = T, xaxt = "n", yaxt = "n", bty = "n",
     col = col_vec2, pch = 16)
axis(side = 1, labels = F)
axis(side = 2, labels = F)
graphics.off()

set.seed(10)
n <- ncol(bm)
idx <- sample(1:n, 5000)
x_vec <- x_vec[idx]
col_vec_subsample <- col_vec[idx]
celltype <- bm$celltype.l2[idx]
ord_idx <- order(x_vec, decreasing = T)
x_vec <- x_vec[ord_idx]
celltype <- celltype[ord_idx]
col_vec_subsample <- col_vec_subsample[ord_idx]

png("../../../out/figures/main/citeseq_bm25_2d_cca_xvec.png", 
    height = 100, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(NA, xlim = range(x_vec), ylim = c(0,1), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:(length(x_vec)-1)){
  polygon(x = c(rep(x_vec[i], 2), rep(x_vec[i+1], 2)),
          y = c(0, 1, 1, 0),
          col = col_vec_subsample[i], border = NA)
}
graphics.off()

set.seed(10)
n <- ncol(bm)
idx <- sample(1:n, 5000)
y_vec <- y_vec[idx]
col_vec_subsample <- col_vec[idx]
celltype <- bm$celltype.l2[idx]
ord_idx <- order(y_vec, decreasing = T)
y_vec <- y_vec[ord_idx]
celltype <- celltype[ord_idx]
col_vec_subsample <- col_vec_subsample[ord_idx]

png("../../../out/figures/main/citeseq_bm25_2d_cca_yvec.png", 
    height = 100, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(NA, xlim = range(y_vec), ylim = c(0,1), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:(length(y_vec)-1)){
  polygon(x = c(rep(y_vec[i], 2), rep(y_vec[i+1], 2)),
          y = c(0, 1, 1, 0),
          col = col_vec_subsample[i], border = NA)
}
graphics.off()

set.seed(10)
n <- ncol(bm)
idx <- sample(1:n, 5000)
projection_vec <- projection_vec[idx,1]
col_vec_subsample <- col_vec[idx]
celltype <- bm$celltype.l2[idx]
ord_idx <- order(projection_vec, decreasing = T)
projection_vec <- projection_vec[ord_idx]
celltype <- celltype[ord_idx]
col_vec_subsample <- col_vec_subsample[ord_idx]

png("../../../out/figures/main/citeseq_bm25_2d_cca_projection.png", 
    height = 100, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(NA, xlim = range(projection_vec), ylim = c(0,1), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:(length(projection_vec)-1)){
  polygon(x = c(rep(projection_vec[i], 2), rep(projection_vec[i+1], 2)),
          y = c(0, 1, 1, 0),
          col = col_vec_subsample[i], border = NA)
}
graphics.off()

############

x_vec <- multiSVD_obj$svd_1$u[,1]
y_vec <- multiSVD_obj$svd_2$u[,1]
x_vec <- 5*(x_vec - stats::median(x_vec))
y_vec <- 5*(y_vec - stats::median(y_vec));
y_vec <- y_vec * sd(x_vec)/sd(y_vec)

pca_res <- stats::prcomp(cbind(x_vec, y_vec))
projection_vec <- cbind(x_vec, y_vec) %*% pca_res$rotation[,1] %*% t(pca_res$rotation[,1])

range(x_vec); range(y_vec)
range(projection_vec[,1]); range(projection_vec[,2])

png("../../../out/figures/main/citeseq_bm25_2d_pca.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(1, 1, 1, 1))
plot(x_vec, y_vec, 
     xaxt = "n", yaxt = "n", bty = "n", asp = T,
     col = col_vec2, pch = 16, xlim = c(-0.04, 0.23), ylim = c(-0.08, 0.12))
points(projection_vec[,1], projection_vec[,2], col = "white", cex = 2, pch = 16)
points(projection_vec[,1], projection_vec[,2], col = col_vec, pch = 16)
axis(side = 1, labels = F)
axis(side = 2, labels = F)
graphics.off()

png("../../../out/figures/main/citeseq_bm25_2d_pca_cleaned.png", 
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(1, 1, 1, 1))
plot(x_vec, y_vec,
     xaxt = "n", yaxt = "n", bty = "n", asp = T,
     col = col_vec2, pch = 16, xlim = c(-0.04, 0.23), ylim = c(-0.08, 0.12))
axis(side = 1, labels = F)
axis(side = 2, labels = F)
graphics.off()

set.seed(10)
n <- ncol(bm)
idx <- sample(1:n, 5000)
x_vec <- x_vec[idx]
col_vec_subsample <- col_vec[idx]
celltype <- bm$celltype.l2[idx]
ord_idx <- order(x_vec, decreasing = T)
x_vec <- x_vec[ord_idx]
celltype <- celltype[ord_idx]
col_vec_subsample <- col_vec_subsample[ord_idx]

png("../../../out/figures/main/citeseq_bm25_2d_pca_xvec.png", 
    height = 100, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(NA, xlim = range(x_vec), ylim = c(0,1), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:(length(x_vec)-1)){
  polygon(x = c(rep(x_vec[i], 2), rep(x_vec[i+1], 2)),
          y = c(0, 1, 1, 0),
          col = col_vec_subsample[i], border = NA)
}
graphics.off()

set.seed(10)
n <- ncol(bm)
idx <- sample(1:n, 5000)
y_vec <- y_vec[idx]
col_vec_subsample <- col_vec[idx]
celltype <- bm$celltype.l2[idx]
ord_idx <- order(y_vec, decreasing = T)
y_vec <- y_vec[ord_idx]
celltype <- celltype[ord_idx]
col_vec_subsample <- col_vec_subsample[ord_idx]

png("../../../out/figures/main/citeseq_bm25_2d_pca_yvec.png", 
    height = 100, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(NA, xlim = range(y_vec), ylim = c(0,1), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:(length(y_vec)-1)){
  polygon(x = c(rep(y_vec[i], 2), rep(y_vec[i+1], 2)),
          y = c(0, 1, 1, 0),
          col = col_vec_subsample[i], border = NA)
}
graphics.off()


set.seed(10)
n <- ncol(bm)
idx <- sample(1:n, 5000)
pca_vec <- pca_res$x[idx,1]
col_vec_subsample <- col_vec[idx]
celltype <- bm$celltype.l2[idx]
ord_idx <- order(pca_vec, decreasing = T)
pca_vec <- pca_vec[ord_idx]
celltype <- celltype[ord_idx]
col_vec_subsample <- col_vec_subsample[ord_idx]

png("../../../out/figures/main/citeseq_bm25_2d_pca_projection.png", 
    height = 100, width = 2500, res = 500, units = "px")
par(mar = c(0.1, 0.1, 0.1, 0.1))
plot(NA, xlim = range(pca_vec), ylim = c(0,1), 
     xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n")
for(i in 1:(length(pca_vec)-1)){
  polygon(x = c(rep(pca_vec[i], 2), rep(pca_vec[i+1], 2)),
          y = c(0, 1, 1, 0),
          col = col_vec_subsample[i], border = NA)
}
graphics.off()