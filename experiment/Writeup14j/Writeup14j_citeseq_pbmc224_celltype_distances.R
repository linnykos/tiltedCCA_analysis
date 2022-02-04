rm(list=ls())
load("../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_preprocessed.RData")
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(pbmc)

# compute low-dimensional embedding
.compute_average_mat <- function(mat, metacell_clustering){
  t(sapply(levels(metacell_clustering), function(celltype){
    idx <- which(metacell_clustering == celltype)
    matrixStats::colSums2(mat[idx,,drop = F])
  }))
}

mat_1 <- pbmc[["pca"]]@cell.embeddings
avgmat_1 <- .compute_average_mat(mat_1, pbmc$celltype.l2)
l2_vec <- apply(avgmat_1, 1, multiomicCCA:::.l2norm)
avgmat_1 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, avgmat_1)

mat_2 <- pbmc[["apca"]]@cell.embeddings
avgmat_2 <- .compute_average_mat(mat_2, pbmc$celltype.l2)
l2_vec <- apply(avgmat_2, 1, multiomicCCA:::.l2norm)
avgmat_2 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, avgmat_2)

distmat_1 <- as.matrix(stats::dist(avgmat_1))
distmat_2 <- as.matrix(stats::dist(avgmat_2))

png("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_distance_rna.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = c(1,6,4,1))
image(multiomicCCA:::.rotate(distmat_1), main = "RNA cosine distances",
      xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n", asp = T)
n <- nrow(distmat_1)
spacing <- 1/(n-1)
y_spacing <- seq(spacing/2, 1, by = 5*spacing)[-1]
x_spacing <- seq(1-spacing/2, 0, by = -5*spacing)[-1]
for(i in y_spacing){
  lines(c(0,1), rep(i,2), lwd = 2, lty = 2)
}
for(i in x_spacing){
  lines(rep(i,2), c(0,1), lwd = 2, lty = 2)
}
axis(side = 2, 
     at = seq(1, 0, by = -spacing), 
     labels = rownames(distmat_1),
     las = 2)
graphics.off()


png("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_distance_adt.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = c(1,6,4,1))
image(multiomicCCA:::.rotate(distmat_2), main = "ADT cosine distances",
      xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n", asp = T)
n <- nrow(distmat_1)
spacing <- 1/(n-1)
y_spacing <- seq(spacing/2, 1, by = 5*spacing)[-1]
x_spacing <- seq(1-spacing/2, 0, by = -5*spacing)[-1]
for(i in y_spacing){
  lines(c(0,1), rep(i,2), lwd = 2, lty = 2)
}
for(i in x_spacing){
  lines(rep(i,2), c(0,1), lwd = 2, lty = 2)
}
axis(side = 2, 
     at = seq(1, 0, by = -spacing), 
     labels = rownames(distmat_1),
     las = 2)
graphics.off()


png("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_distance_rna_logged.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = c(1,6,4,1))
tmp <- exp(distmat_1)
image(multiomicCCA:::.rotate(tmp), main = "RNA cosine distances\n(Exponentiated)",
      xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n", asp = T)
n <- nrow(distmat_1)
spacing <- 1/(n-1)
y_spacing <- seq(spacing/2, 1, by = 5*spacing)[-1]
x_spacing <- seq(1-spacing/2, 0, by = -5*spacing)[-1]
for(i in y_spacing){
  lines(c(0,1), rep(i,2), lwd = 2, lty = 2)
}
for(i in x_spacing){
  lines(rep(i,2), c(0,1), lwd = 2, lty = 2)
}
axis(side = 2, 
     at = seq(1, 0, by = -spacing), 
     labels = rownames(distmat_1),
     las = 2)
graphics.off()


png("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_distance_adt_logged.png",
    height = 2500, width = 2500, units = "px", res = 300)
par(mar = c(1,6,4,1))
tmp <- exp(distmat_2)
image(multiomicCCA:::.rotate(tmp), main = "ADT cosine distances\n(Exponentiated)",
      xlab = "", ylab = "", bty = "n", xaxt = "n", yaxt = "n", asp = T)
n <- nrow(distmat_1)
spacing <- 1/(n-1)
y_spacing <- seq(spacing/2, 1, by = 5*spacing)[-1]
x_spacing <- seq(1-spacing/2, 0, by = -5*spacing)[-1]
for(i in y_spacing){
  lines(c(0,1), rep(i,2), lwd = 2, lty = 2)
}
for(i in x_spacing){
  lines(rep(i,2), c(0,1), lwd = 2, lty = 2)
}
axis(side = 2, 
     at = seq(1, 0, by = -spacing), 
     labels = rownames(distmat_1),
     las = 2)
graphics.off()