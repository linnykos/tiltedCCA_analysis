rm(list=ls())
set.seed(10)

library(Seurat); library(Signac)
load("../../out/Writeup9_10x_pbmc_dimred_all.RData")

head(pbmc@meta.data)
pbmc <- Seurat::FindMultiModalNeighbors(pbmc, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
pbmc <- Seurat::RunUMAP(pbmc, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc <- Seurat::FindClusters(pbmc, graph.name = "wsnn", algorithm = 3, verbose = FALSE)
pbmc$celltype <- Seurat::Idents(pbmc)

mat <- pbmc[["SCT"]]@scale.data
set.seed(10)
kmeans_res <- stats::kmeans(t(mat), centers = 10)
sample_idx <- unlist(lapply(1:10, function(x){
  idx <- which(kmeans_res$cluster == x)
  if(length(idx) < 200){return(idx)}
  sample(idx, 200)
}))
table(pbmc$celltype[sample_idx])
mat <- mat[,sample_idx]

K <- 15
set.seed(10)
kmeans_res <- stats::kmeans(mat, centers = K)
var_order <- as.numeric(unlist(lapply(1:K, function(x){
  idx <- which(kmeans_res$cluster == x)
  sample(idx, size = round(length(idx)*2/3))
})))

tmp <- lapply(unique(pbmc$celltype), function(x){
  sample_idx[which(pbmc$celltype[sample_idx] == x)]
})
type_length <- sapply(tmp, length)
sample_idx2 <- unlist(tmp)

mat1 <- t(as.matrix(pbmc[["RNA"]]@counts)[var_order, sample_idx2])
mat2 <- t(as.matrix(pbmc[["SCT"]]@data)[var_order, sample_idx2])
mat3 <- t(pbmc[["SCT"]]@scale.data[var_order, sample_idx2])

png("../../out/Writeup9_10x_pbmc_rna1.png", height = 2500, width = 2500, units = "px", res = 300)
plot_heatmat(mat1, xlab = "Genes", ylab = "Cells", main = "PBMC (10x): RNA Raw counts")
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

png("../../out/Writeup9_10x_pbmc_rna2.png", height = 2500, width = 2500, units = "px", res = 300)
plot_heatmat(mat2, xlab = "Genes", ylab = "Cells", main = "PBMC (10x): RNA SCT-norm. counts")
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

png("../../out/Writeup9_10x_pbmc_rna3.png", height = 2500, width = 2500, units = "px", res = 300)
plot_heatmat(mat3, reserve_zero = F, luminosity = F, xlab = "Genes", ylab = "Cells", main = "PBMC (10x): RNA Scaled counts")
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

quantile_vec1 <- stats::quantile(mat1, probs = seq(0,1,length.out = 21))
quantile_vec2 <- stats::quantile(mat2, probs = seq(0,1,length.out = 21))
quantile_vec3 <- stats::quantile(mat3, probs = seq(0,1,length.out = 21))
quantile_vec1; quantile_vec2; quantile_vec3
length(which(mat1 == 0))/prod(dim(mat1))

png("../../out/Writeup9_10x_pbmc_rna3_hist.png", height = 1500, width = 2500, units = "px", res = 300)
graphics::hist(mat3, breaks = 50, col = "gray", 
               xlab = "Value",
               main = "PBMC (10x): RNA Scaled counts")
graphics.off()

##############################################

mat <- pbmc[["ATAC"]]@scale.data[,sample_idx]
row_vec <- rowMeans(mat)
var_idx <- order(abs(row_vec), decreasing = T)[1:5000]
K <- 15
set.seed(10)
kmeans_res <- stats::kmeans(mat[var_idx,], centers = K)
var_order <- var_idx[as.numeric(unlist(lapply(1:K, function(x){
  idx <- which(kmeans_res$cluster == x)
  sample(idx, size = round(length(idx)*2000/length(var_idx)))
})))]

mat1 <- t(as.matrix(pbmc[["ATAC"]]@counts)[var_order, sample_idx2])
mat2 <- t(as.matrix(pbmc[["ATAC"]]@data)[var_order, sample_idx2])
mat3 <- t(pbmc[["ATAC"]]@scale.data[var_order, sample_idx2])

png("../../out/Writeup9_10x_pbmc_atac1.png", height = 2500, width = 2500, units = "px", res = 300)
plot_heatmat(mat1, xlab = "ATAC", ylab = "Cells", main = "PBMC (10x): ATAC Raw counts")
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

png("../../out/Writeup9_10x_pbmc_atac2.png", height = 2500, width = 2500, units = "px", res = 300)
plot_heatmat(mat2, xlab = "ATAC", ylab = "Cells", main = "PBMC (10x): ATAC TF-IDF counts")
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

png("../../out/Writeup9_10x_pbmc_atac3.png", height = 2500, width = 2500, units = "px", res = 300)
plot_heatmat(mat3, reserve_zero = F, luminosity = F, xlab = "ATAC", ylab = "Cells", main = "PBMC (10x): ATAC Scaled counts")
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

quantile_vec1 <- stats::quantile(mat1, probs = seq(0,1,length.out = 21))
quantile_vec2 <- stats::quantile(mat2, probs = seq(0,1,length.out = 21))
quantile_vec3 <- stats::quantile(mat3, probs = seq(0,1,length.out = 21))
quantile_vec1; quantile_vec2; quantile_vec3
length(which(mat1 == 0))/prod(dim(mat1))
length(pbmc[["ATAC"]]@counts@i)/prod(dim(pbmc[["ATAC"]]@counts))
length(which(pbmc[["ATAC"]]@counts@x==1))/length(pbmc[["ATAC"]]@counts@x)
length(which(pbmc[["ATAC"]]@counts@x==2))/length(pbmc[["ATAC"]]@counts@x)
length(which(pbmc[["ATAC"]]@counts@x==3))/length(pbmc[["ATAC"]]@counts@x)
length(which(pbmc[["ATAC"]]@counts@x==4))/length(pbmc[["ATAC"]]@counts@x)
length(which(pbmc[["ATAC"]]@counts@x<=4))/length(pbmc[["ATAC"]]@counts@x)


png("../../out/Writeup9_10x_pbmc_atac3_hist.png", height = 1500, width = 2500, units = "px", res = 300)
graphics::hist(mat3, breaks = 50, col = "gray", 
               xlab = "Value",
               main = "PBMC (10x): ATAC Scaled counts")
graphics.off()

