rm(list=ls()); set.seed(10); # gcinfo(TRUE)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

date_of_run <- Sys.time(); session_info <- sessionInfo()
load("../../out/Writeup10_10x_pbmc_preprocess4.RData"); pbmc[["ATAC"]] = NULL
load("../../out/Writeup10_10x_pbmc_dcca_metacells2.RData")

#############

res <- dcca_decomposition(dcca_res, rank_c = K)
rownames(res$common_mat_1) <- colnames(pbmc)

1-sqrt((1-dcca_res$cca_obj)/(1+dcca_res$cca_obj))

#######################################
.l2norm <- function(x){sqrt(sum(x^2))}
rank_1 <- ncol(res$distinct_score_1); rank_2 <- ncol(res$distinct_score_2)
tmp1 <- sapply(1:ncol(res$common_score), function(j){
  .l2norm(res$common_score[,j])/(.l2norm(res$common_score[,j]) + .l2norm(res$distinct_score_1[,j]))
})

tmp4 <- multiomicCCA:::.svd_truncated(res$common_mat_1 + res$distinct_mat_1, rank_1)
tmp5 <- multiomicCCA:::.svd_truncated(res$common_mat_2 + res$distinct_mat_2, rank_2)

png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_values.png", height = 800, width = 2000, units = "px", res = 300)
par(mfrow = c(1,3))
# cca values
plot(dcca_res$cca_obj, main = "CCA objective", xlab = "Canonical dimension", 
     ylab = "Canonical correlation", pch = 16)

# amount of variability, common vs distincts
plot(NA, xlim = c(1,length(tmp1)), ylim = c(0,1), main = "Ratio of canonical\nscores' signal",
     xlab = "Canonical dimension", ylab = "Ratio of signals (common/total)")
points(tmp1, pch = 16)
for(i in 1:length(tmp1)){
  lines(rep(i,2), c(0,tmp1[i]), col = 1)
  lines(rep(i,2), c(1,tmp1[i]), col = 2)
}

# amount of variability, rna vs atac
ylim <- range(c(tmp4$d, tmp5$d))
plot(NA, xlim = c(1,max(c(rank_1,rank_2))), ylim = ylim, main = "Singular value of modalities",
     xlab = "Singular dimension", ylab = "Singular value")
points(tmp4$d, pch = 16)
points(tmp5$d, col = 2, pch = 16)
legend("topright", c("RNA", "ATAC"), bty="n", fill=c(1,2))
graphics.off()

###########################################

# try computing RNA-vs-ATAC weights
# for each cell, compute the proportion in distinct vs common
n <- nrow(res$common_mat_1)
.l2norm <- function(x){sqrt(sum(x^2))}
weight_mat <- t(sapply(1:n, function(i){
  if(i %% floor(n/10) == 0) cat('*')
  vec <- c(.l2norm(res$distinct_mat_1[i,])/(.l2norm(res$common_mat_1[i,])+.l2norm(res$distinct_mat_1[i,])),
           .l2norm(res$distinct_mat_2[i,])/(.l2norm(res$common_mat_2[i,])+.l2norm(res$distinct_mat_2[i,])))
  names(vec) <- c("RNA", "ATAC")
  vec
}))

quantile(weight_mat[,1]); quantile(weight_mat[,2])
weight_mat <- data.frame(RNA = weight_mat[,1], ATAC = weight_mat[,2], cell_type = pbmc@meta.data[,"predicted.id"])

# compute group-wise weight matrix
uniq_types <- unique(weight_mat$cell_type)
summary_mat <- t(sapply(uniq_types, function(i){
  idx <- which(weight_mat$cell_type == i)
  colMeans(weight_mat[idx,1:2,drop = F])
}))
summary_mat <- data.frame(RNA = summary_mat[,1], ATAC = summary_mat[,2], cell_type = uniq_types)

# order the weight matrix
diff_val <- summary_mat[,1] - summary_mat[,2]
summary_mat <- summary_mat[order(diff_val, decreasing = T),]

# take the colors
col_vec <- scales::hue_pal()(length(uniq_types))
type_col_mat <- cbind(as.character(sort(uniq_types)), col_vec)

# draw the plot
len <- length(uniq_types)
png("../../out/figures/Writeup10/Writeup10_pbmc_10x_dcca_weights.png", height = 1500, width = 2500, units = "px", res = 300)
par(mar = c(6, 4, 3.5, 0.5))
plot(NA, xlim = c(0.5, 1.5*(len-1)), ylim = c(-1,1), xlab = "", ylab = "Ratio explained",
     xaxt = "n", main = "RNA distinct vs. Common vs. ATAC distinct\nexplained variability")
title(xlab = "Cell type", mgp = c(4,1,0))

y_vec <- c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75)
for(i in 1:length(y_vec)){
  lines(c(-10*len, 10*len), rep(y_vec[i], 2), col = "red", lwd = 0.75, lty = 3)
}

for(i in 1:len){
  xlim <- c((i-1)*1.5, (i-1)*1.5+1)
  
  set.seed(10)
  idx <- which(weight_mat$cell_type == summary_mat[i,3])
  if(length(idx) > 50){ idx <- idx[sample(1:length(idx), 50)] }
  
  col <- type_col_mat[which(type_col_mat[,1] == summary_mat[i,3]), 2]
  graphics::rect(xlim[1], -(1-summary_mat[i,2]), xlim[2], 1-summary_mat[i,1], col = col)
  offsets <- seq(0, 1, length.out = length(idx)+10)
  offsets <- offsets[-c(1:5, (length(offsets)-4):(length(offsets)))]
  
  for(j in 1:length(idx)){
    lines(rep(xlim[1]+offsets[j], 2), c(1-weight_mat[idx[j],1], -(1-weight_mat[idx[j],2])),
          col = "gray", lwd = 0.3)
    points(rep(xlim[1]+offsets[j], 2), c(1-weight_mat[idx[j],1], -(1-weight_mat[idx[j],2])),
           col = "gray", cex = 0.2, pch = 16)
  }
}
lines(c(-2*len, 2*len), rep(0,2), col = "black", lwd = 3, lty = 2)

axis(side = 1, at = (((1:len)-1)*1.5+1.5/2), labels = FALSE)
text(x = (((1:len)-1)*1.5+1.5/2),
     y = par("usr")[3]-0.1,
     labels = as.character(summary_mat[,3]),
     xpd = NA, srt = 45, adj = 0.965, cex = 0.8)

graphics.off()

# var_list <- ls()
# var_list <- var_list[!var_list %in% "weight_mat"]
# rm(list = var_list)
# 
# save.image("../../out/Writeup10_10x_pbmc_weight.RData")

