rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/Writeup14p/Writeup14p_10x_greenleaf_tcca.RData")

rna_mat <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_mat <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_mat <- rna_mat[cell_idx[order_vec],]
atac_mat <- atac_mat[cell_idx[order_vec],]

save(rna_mat, atac_mat,
     file = "../../../../out/Writeup14p/tmp.RData")

#########################

load("../../out/experiment/Writeup14p/tmp.RData")

mentioned_genes <- c("ALDH2", "APOE", "AQP4", "ASCL1", "BHLHE22",
                     "BHLHE40", "C16orf89", "CAV2", "DLX2", "DOK5",
                     "DUSP1", "EOMES", "ETV4", "FOS", "FOXJ1", "GLI3",
                     "HAS2", "HES1", "HES4", "HSPA1A", "HSPA1B",
                     "ID3", "IGFBP7", "JUN", "KIF1A", "LIMCH1",
                     "MBP", "MEF2C", "NEUROD1", "NEUROD2", "NEUROD4",
                     "NEUROD6", "NEUROG1", "NEUROG2", "NFIA", "NFIB", "NFIC",
                     "NHLH1", "NR2F1", "PAX6", "RFX4", "RUNX1",
                     "OLIG1", "OLIG2", "SOX2", "SOX3", "SOX6",
                     "SOX9", "SOX10", "SOX21", "SPARCL1", "SNCB", "TBX",
                     "TNC", "TOP2A", "TRB1", "WNT11",
                     "SATB2")
mentioned_genes <- sort(intersect(mentioned_genes, colnames(rna_mat)))

rna_mat_subset <- rna_mat[,mentioned_genes]
atac_mat_subset <- atac_mat[,mentioned_genes]
plot(rna_mat_subset[,"EOMES"], pch = 16)
points(atac_mat_subset[,"EOMES"], pch = 16, col = 2)

rna_mat2 <- t(scale(t(rna_mat)))
atac_mat2 <- t(scale(t(atac_mat)))
rna_mat_subset2 <- rna_mat2[,mentioned_genes]
atac_mat_subset2 <- atac_mat2[,mentioned_genes]
plot(rna_mat_subset2[,"EOMES"], pch = 16)
points(atac_mat_subset2[,"EOMES"], pch = 16, col = 2)

diff_mat <- (rna_mat_subset2 - atac_mat_subset2)/abs(rna_mat_subset2)
diff_matb <- pmax(pmin(diff_mat, 15), -15)
plot(diff_matb[,"EOMES"])
plot(diff_matb[,"EOMES"], ylim = c(-1,1)); lines(c(0,1e5), rep(0,2), col = 2, lwd = 2, lty = 2)
diff_mat2 <- abs(diff_mat)
diff_mat2 <- pmin(diff_mat2, 15)
plot(diff_mat2[,"EOMES"])

n <- nrow(diff_matb)
vec <- diff_matb[,"EOMES"]
bw <- 500
smooth_vec <- sapply(1:n, function(i){
  mean(vec[max(1,i-bw):min(n,i+bw)])
})
plot(smooth_vec)

tmp_df <- data.frame(y = smooth_vec, x = 1:n)
reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
x_vec <- 1:n
y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
smooth_vec2 <- y_vec$Estimation[,"Pred"]
plot(smooth_vec2)

rna_vec <- rna_mat_subset2[,"EOMES"]
tmp_df <- data.frame(y = rna_vec, x = 1:n)
reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
x_vec <- 1:n
y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
smooth_rna_vec <- y_vec$Estimation[,"Pred"]
plot(smooth_rna_vec)
atac_vec <- smooth_rna_vec - abs(smooth_rna_vec)*smooth_vec2

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))
plot(x = atac_vec, y = smooth_rna_vec, col = color_vec, pch = 16, asp = T)

#####################################

diff_mat <- rna_mat_subset2 - atac_mat_subset2
plot(diff_mat[,"EOMES"]); lines(c(0,1e5), rep(0,2), col = 2, lwd = 2, lty = 2)
diff_vec <- diff_mat[,"EOMES"]

tmp_df <- data.frame(y = diff_vec, x = 1:n)
reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
x_vec <- 1:n
y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
smooth_diff_vec <- y_vec$Estimation[,"Pred"]
plot(smooth_diff_vec)

atac_vec <- smooth_rna_vec - smooth_diff_vec
plot(smooth_rna_vec)
points(atac_vec, col = 2)

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))
plot(x = atac_vec, y = smooth_rna_vec, col = color_vec, pch = 16, asp = T)
plot(x = smooth_diff_vec, y = smooth_rna_vec, col = color_vec, pch = 16)

########################
########################
########################
########################

diff_mat <- rna_mat_subset2 - atac_mat_subset2
diff_vec <- diff_mat[,"ASCL1"]
plot(diff_vec); lines(c(0,1e5), rep(0,2), col = 2, lwd = 2, lty = 2)

bw <- 1000
smooth_vec <- sapply(1:n, function(i){
  mean(diff_vec[max(1,i-bw):min(n,i+bw)])
})
plot(smooth_vec)

tmp_df <- data.frame(y = smooth_vec, x = 1:n)
reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
x_vec <- 1:n
y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
smooth_diff_vec <- y_vec$Estimation[,"Pred"]
plot(smooth_diff_vec)

rna_vec <- rna_mat_subset2[,"ASCL1"]
tmp_df <- data.frame(y = rna_vec, x = 1:n)
reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
x_vec <- 1:n
y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
smooth_rna_vec <- y_vec$Estimation[,"Pred"]
plot(smooth_rna_vec)

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))
plot(x = smooth_diff_vec, y = smooth_rna_vec, col = color_vec, pch = 16)


##################################

suggested_pipeline <- function(rna_mat_subset2,
                               atac_mat_subset2,
                               gene,
                               bw = 1000){
  # diff_mat <- (rna_mat_subset2 - atac_mat_subset2)/abs(rna_mat_subset2)
  # try this for now cuz i'm crazy
  diff_mat <- abs(rna_mat_subset2 - atac_mat_subset2)
  diff_vec <- diff_mat[,gene]
  diff_vec <- as.numeric(scale(diff_vec))
  n <- length(diff_vec)
  plot(diff_vec, main = "diff_vec"); lines(c(-1e6,1e6), rep(0,2), col = 2, lwd = 2, lty = 2)
  
  bw <- 1000
  smooth_vec <- sapply(1:n, function(i){
    mean(diff_vec[max(1,i-bw):min(n,i+bw)])
  })
  plot(smooth_vec, main = "smooth_vec"); lines(c(-1e6,1e6), rep(0,2), col = 2, lwd = 2, lty = 2)
  
  tmp_df <- data.frame(y = smooth_vec, x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  smooth_diff_vec <- y_vec$Estimation[,"Pred"]
  plot(smooth_diff_vec, main = "smooth_diff_vec"); lines(c(-1e6,1e6), rep(0,2), col = 2, lwd = 2, lty = 2)
  
  rna_vec <- rna_mat_subset2[,gene]
  bw <- 1000
  smooth_rna_vec <- sapply(1:n, function(i){
    mean(rna_vec[max(1,i-bw):min(n,i+bw)])
  })
  plot(smooth_rna_vec, main = "smooth_rna_vec"); lines(c(-1e6,1e6), rep(0,2), col = 2, lwd = 2, lty = 2)
  
  tmp_df <- data.frame(y = smooth_rna_vec, x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  smooth_rna_vec2 <- y_vec$Estimation[,"Pred"]
  plot(smooth_rna_vec2, main = "smooth_rna_vec2"); lines(c(-1e6,1e6), rep(0,2), col = 2, lwd = 2, lty = 2)
  
  base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
  color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))
  plot(x = smooth_diff_vec, y = smooth_rna_vec2, col = color_vec, pch = 16, main = "leafplot")
}

suggested_pipeline(rna_mat_subset2, atac_mat_subset2, gene = "EOMES")
suggested_pipeline(rna_mat_subset2, atac_mat_subset2, gene = "SOX2")
suggested_pipeline(rna_mat_subset2, atac_mat_subset2, gene = "GLI3")
suggested_pipeline(rna_mat_subset2, atac_mat_subset2, gene = "WNT11")
suggested_pipeline(rna_mat_subset2, atac_mat_subset2, gene = "SATB2")
suggested_pipeline(rna_mat_subset2, atac_mat_subset2, gene = "SOX21")

suggested_pipeline(rna_mat_subset2, atac_mat_subset2, gene = "NEUROD1")
