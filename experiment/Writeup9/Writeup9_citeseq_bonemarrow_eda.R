rm(list=ls())
set.seed(10)

# from https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html
date_of_run <- Sys.time()
session_info <- sessionInfo()

library(Seurat)
load("../../out/Writeup9_citeseq_bonemarrow_dimred_all.RData")

# let's do some fast and loose subsampling to get the data more tractable
# bm <- Seurat::FindMultiModalNeighbors(
#   bm, reduction.list = list("pca", "apca"), 
#   dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight"
# )
# 
# bm <- Seurat::RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# bm <- Seurat::FindClusters(bm, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = FALSE)

head(bm@meta.data)
length(unique(bm@meta.data[,"celltype.l1"])); length(unique(bm@meta.data[,"celltype.l2"]))
# length(unique(bm@meta.data[,"seurat_clusters"]))
table(bm@meta.data[,"celltype.l2"])

mat <- bm[["RNA"]]@scale.data
set.seed(10)
kmeans_res <- stats::kmeans(t(mat), centers = 10)
sample_idx <- unlist(lapply(1:10, function(x){
  idx <- which(kmeans_res$cluster == x)
  if(length(idx) < 200){return(idx)}
  sample(idx, 200)
}))
table(bm@meta.data[sample_idx,"celltype.l2"])
mat <- mat[,sample_idx]

K <- 15
set.seed(10)
kmeans_res <- stats::kmeans(mat, centers = K)
var_order <- unlist(lapply(1:K, function(x){
  which(kmeans_res$cluster == x)
}))

tmp <- lapply(unique(bm@meta.data[,"celltype.l2"]), function(x){
  sample_idx[which(bm@meta.data[sample_idx,"celltype.l2"] == x)]
})
type_length <- sapply(tmp, length)
sample_idx2 <- unlist(tmp)

mat1 <- t(as.matrix(bm[["RNA"]]@counts)[var_order, sample_idx2])
mat2 <- t(as.matrix(bm[["RNA"]]@data)[var_order, sample_idx2])
mat3 <- t(bm[["RNA"]]@scale.data[var_order, sample_idx2])

png("../../out/Writeup9_citeseq_bonemarrow_rna1.png", height = 2500, width = 2500, units = "px", res = 300)
plot_heatmat(mat1, xlab = "Genes", ylab = "Cells", main = "Bone marrow (CITE-seq): RNA Raw counts")
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

png("../../out/Writeup9_citeseq_bonemarrow_rna2.png", height = 2500, width = 2500, units = "px", res = 300)
plot_heatmat(mat2, xlab = "Genes", ylab = "Cells", main = "Bone marrow (CITE-seq): RNA Log-normalized counts")
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

png("../../out/Writeup9_citeseq_bonemarrow_rna3.png", height = 2500, width = 2500, units = "px", res = 300)
plot_heatmat(mat3, reserve_zero = F, luminosity = F, xlab = "Genes", ylab = "Cells", main = "Bone marrow (CITE-seq): RNA Scaled counts")
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

png("../../out/Writeup9_citeseq_bonemarrow_rna3_hist.png", height = 1500, width = 2500, units = "px", res = 300)
graphics::hist(mat3, breaks = 50, col = "gray", 
               xlab = "Value",
               main = "Bone marrow (CITE-seq): RNA Scaled counts")
graphics.off()

################################################


mat3 <- t(bm[["ADT"]]@scale.data[, sample_idx2])
var_order <- stats::hclust(stats::dist(t(mat3)))$order

mat1 <- t(as.matrix(bm[["ADT"]]@counts)[var_order, sample_idx2])
mat2 <- t(as.matrix(bm[["ADT"]]@data)[var_order, sample_idx2])
mat3 <- t(bm[["ADT"]]@scale.data[var_order, sample_idx2])

png("../../out/Writeup9_citeseq_bonemarrow_adt1.png", height = 2500, width = 1500, units = "px", res = 300)
plot_heatmat(mat1, xlab = "Antibody", ylab = "Cells", main = "Bone marrow (CITE-seq): ADT Raw counts",
             asp = 25/10)
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

png("../../out/Writeup9_citeseq_bonemarrow_adt2.png", height = 2500, width = 1500, units = "px", res = 300)
plot_heatmat(mat2, xlab = "Antibody", ylab = "Cells", main = "Bone marrow (CITE-seq): ADT CLR-norm. counts",
             asp = 25/10)
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

png("../../out/Writeup9_citeseq_bonemarrow_adt3.png", height = 2500, width = 1500, units = "px", res = 300)
plot_heatmat(mat3, reserve_zero = F, luminosity = F, xlab = "Antibody", ylab = "Cells", main = "Bone marrow (CITE-seq): ADT Scaled counts",
             asp = 25/10)
for(i in 1:(length(type_length)-1)){
  n <- nrow(mat1)
  y_val <- (n-sum(type_length[1:i]))/n
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, col = "white")
  graphics::lines(c(0,1), rep(y_val, 2), lwd = 2, lty = 2)
}
graphics.off()

png("../../out/Writeup9_citeseq_bonemarrow_adt3_hist.png", height = 1500, width = 2500, units = "px", res = 300)
graphics::hist(mat3, breaks = 50, col = "gray", 
               xlab = "Value",
               main = "Bone marrow (CITE-seq): ADT Scaled counts")
graphics.off()

quantile_vec1 <- stats::quantile(mat1, probs = seq(0,1,length.out = 21))
quantile_vec2 <- stats::quantile(mat2, probs = seq(0,1,length.out = 21))
quantile_vec3 <- stats::quantile(mat3, probs = seq(0,1,length.out = 21))
quantile_vec1; quantile_vec2; quantile_vec3
length(which(mat1 == 0))/prod(dim(mat1))
