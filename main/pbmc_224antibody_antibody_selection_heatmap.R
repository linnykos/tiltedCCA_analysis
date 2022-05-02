rm(list=ls())
load("../../../out/main/citeseq_pbmc224_varSelect.RData")
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")
source("pbmc_224antibody_colorPalette.R")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

celltype_vec <- factor(pbmc$celltype.l2)
celltype_list <- list("CD4" = c("CD4 CTL", "CD4 Naive", "CD4 Proliferating", "CD4 TCM", "CD4 TEM"),
                      "CD8" = c("CD8 Naive", "CD8 Proliferating", "CD8 TCM", "CD8 TEM"),
                      "OtherT" = c("dnT", "gdT", "MAIT"),
                      "NK" = c("NK", "NK Proliferating", "NK_CD56bright"),
                      "Monocyte" = c("CD14 Mono", "CD16 Mono"),
                      "Dendritic" = c("ASDC", "cDC1", "cDC2", "pDC"),
                      "B" = c("B intermediate", "B memory", "B naive"),
                      "Other" = c("Eryth", "HSPC", "ILC", "Platelet"))
celltype_all <- unlist(celltype_list)
names(celltype_all) <- NULL

protein_mat <- multiSVD_obj$common_mat_2 + multiSVD_obj$distinct_mat_2
protein_mat <- protein_mat[,variable_selection_res$selected_variables]
protein_avg <- t(sapply(celltype_all, function(celltype){
  idx <- which(celltype_vec == celltype)
  colMeans(protein_mat[idx,,drop = F])
}))
protein_avg <- scale(protein_avg)

hclust_res <- stats::hclust(stats::dist(t(protein_avg)))
protein_avg2 <- protein_avg[,hclust_res$order]

break_vec <- c(seq(min(protein_avg2), 0, length.out = 11), seq(0, max(protein_avg2), length.out = 11)[-1])
break_vec[1] <- break_vec[1]-1
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+1
col_vec1 <- grDevices::colorRampPalette(c(rgb(184, 54, 220, maxColorValue = 255), "white"))(11)[-1]
col_vec2 <- grDevices::colorRampPalette(c("white",  rgb(235, 134, 47, maxColorValue = 255)))(11)[-1]
col_vec <- c(col_vec1, col_vec2)

png("../../../out/figures/main/citeseq_pbmc224_varSelect_protein-everything_heatmap.png",
    height = 3500, width = 3500*ncol(protein_avg2)/nrow(protein_avg2), res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(protein_avg2), 
      asp = (nrow(protein_avg2)-1)/(ncol(protein_avg2)-1),
      main = "", xlab = "", ylab = "",
      xaxt = "n", yaxt = "n", bty = "n",
      breaks = break_vec, col = col_vec)

# draw lines
n <- nrow(protein_avg2)
halfspacing <- 1/(2*(n-1))
num_types <- sapply(celltype_list, length)
for(i in 1:(length(num_types)-1)){
  y_val <- 1-((2*halfspacing)*(sum(num_types[1:i])-1)+halfspacing)
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 5, col = "white")
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 4, lty = 2, col = "black")
}
graphics.off()

####################

protein_distinct_mat <- multiSVD_obj$distinct_mat_2
protein_distinct_mat <- protein_distinct_mat[,variable_selection_res$selected_variables]
protein_distinct_avg <- t(sapply(celltype_all, function(celltype){
  idx <- which(celltype_vec == celltype)
  colMeans(protein_distinct_mat[idx,,drop = F])
}))
protein_distinct_avg <- scale(protein_distinct_avg)
protein_distinct_avg2 <- protein_distinct_avg[,hclust_res$order]

break_vec <- c(seq(min(protein_distinct_avg2), 0, length.out = 11), seq(0, max(protein_distinct_avg2), length.out = 11)[-1])
break_vec[1] <- break_vec[1]-1
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+1
col_vec1 <- grDevices::colorRampPalette(c(rgb(184, 54, 220, maxColorValue = 255), "white"))(11)[-1]
col_vec2 <- grDevices::colorRampPalette(c("white",  rgb(235, 134, 47, maxColorValue = 255)))(11)[-1]
col_vec <- c(col_vec1, col_vec2)

png("../../../out/figures/main/citeseq_pbmc224_varSelect_protein-distinct_heatmap.png",
    height = 3500, width = 3500*ncol(protein_distinct_avg2)/nrow(protein_distinct_avg2), res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(protein_distinct_avg2), 
      asp = (nrow(protein_distinct_avg2)-1)/(ncol(protein_distinct_avg2)-1),
      main = "", xlab = "", ylab = "",
      xaxt = "n", yaxt = "n", bty = "n",
      breaks = break_vec, col = col_vec)
graphics.off()
