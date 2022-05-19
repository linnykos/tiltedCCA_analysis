rm(list=ls())
load("../../../out/main/citeseq_bm25_tcca.RData")

library(Seurat)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
adt_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = rna_common[i,],
                   adt = adt_pred[i,])
  lm_res <- stats::lm(rna ~ adt, data = df)
  summary(lm_res)$r.squared
})

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec^val, rank(alignment_vec))
})
bm$alignment <- alignment_vec^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(bm$alignment), max(bm$alignment), length.out = num_color)
color_vec <- sapply(bm$alignment, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

png(paste0("../../../out/figures/main/citeseq_bm25_tcca_steadystate-full.png"),
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = c(4,4,4,0.5))
plot(x = bm[["common_tcca"]]@cell.embeddings[,1],
     y = bm[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     xlab = colnames(bm[["common_tcca"]]@cell.embeddings)[1],
     ylab = colnames(bm[["common_tcca"]]@cell.embeddings)[2],
     main = paste0("Human BM (CITE-Seq, RNA+ADT)\nAlignment between ADT and common RNA"),
     xaxt = "n", yaxt = "n", bty = "n")
axis(side = 1)
axis(side = 2)
graphics.off()

###########

Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@data[Seurat::VariableFeatures(object = bm),])
mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}
alignment_vec_alt <- sapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  df <- data.frame(rna = mat_1b[i,],
                   adt = adt_pred[i,])
  lm_res <- stats::lm(rna ~ adt, data = df)
  summary(lm_res)$r.squared
})

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec_alt^val, rank(alignment_vec_alt))
})
bm$alignment_alt <- alignment_vec_alt^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(bm$alignment_alt), max(bm$alignment_alt), length.out = num_color)
color_vec <- sapply(bm$alignment_alt, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

png(paste0("../../../out/figures/main/citeseq_bm25_tcca_steadystate_alt-full.png"),
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = c(4,4,4,0.5))
plot(x = bm[["common_tcca"]]@cell.embeddings[,1],
     y = bm[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     xlab = colnames(bm[["common_tcca"]]@cell.embeddings)[1],
     ylab = colnames(bm[["common_tcca"]]@cell.embeddings)[2],
     main = paste0("Human BM (CITE-Seq, RNA+ADT)\nAlignment between ADT and RNA"),
     xaxt = "n", yaxt = "n", bty = "n")
axis(side = 1)
axis(side = 2)
graphics.off()