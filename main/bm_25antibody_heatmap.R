rm(list=ls())
library(Seurat)
load("../../../out/main/citeseq_bm25_tcca.RData")
source("bm_25antibody_colorPalette.R")

celltype_vec <- factor(bm$celltype.l2)
celltype_list <- list("CD4" = c("CD4 Memory", "CD4 Naive", "Treg"),
                      "CD8" = c("CD8 Effector_1", "CD8 Effector_2", "CD8 Memory_1", "CD8 Memory_2", "CD8 Naive", "gdT", "MAIT"), 
                      "NK cells" = c("CD56 bright NK", "NK"),
                      "B" = c("Memory B", "Naive B", "Plasmablast"),
                      "Myeloid" = c("CD14 Mono", "CD16 Mono", "cDC2", "pDC"),
                      "Progenitor" = c("GMP", "HSC", "LMPP", "Prog_DC", "Prog_Mk", "Prog_RBC", "Prog_B 1", "Prog_B 2"))
celltype_all <- unlist(celltype_list)
names(celltype_all) <- NULL

.processing_func <- function(mat, val = 0.05){
  mat <- scale(mat, center = T, scale = F)
  pmin(pmax(mat, quantile(mat, probs = val)), quantile(mat,probs = 1-val))
}

rna_mat <- multiSVD_obj$svd_1$u  %*% diag(multiSVD_obj$svd_1$d)
rna_mat <- .processing_func(rna_mat)
rna_avg <- t(sapply(celltype_all, function(celltype){
  idx <- which(celltype_vec == celltype)
  colMeans(rna_mat[idx,,drop = F])
}))
rna_avg <- rna_avg[,1:10]

break_vec <- c(seq(min(rna_avg), 0, length.out = 11), seq(0, max(rna_avg), length.out = 11)[-1])
break_vec[1] <- break_vec[1]-1
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+1
base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
col_vec <- rev(grDevices::colorRampPalette(base_palette)(length(break_vec)-1))

png("../../../out/figures/main/bm_25antibody_rna-heatmap.png",
    height = 3500, width = 3500*ncol(rna_avg)/nrow(rna_avg), res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(rna_avg), 
      asp = (nrow(rna_avg)-1)/(ncol(rna_avg)-1),
      main = "", xlab = "", ylab = "",
      xaxt = "n", yaxt = "n", bty = "n",
      breaks = break_vec, col = col_vec)

# draw lines
n <- nrow(rna_avg)
halfspacing <- 1/(2*(n-1))
num_types <- sapply(celltype_list, length)
for(i in 1:(length(num_types)-1)){
  y_val <- 1-((2*halfspacing)*(sum(num_types[1:i])-1)+halfspacing)
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 5, col = "white")
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 4, lty = 2, col = "black")
}
graphics.off()

##############

adt_mat <- multiSVD_obj$svd_2$u %*% diag(multiSVD_obj$svd_2$d)
adt_mat <- .processing_func(adt_mat)
adt_avg <- t(sapply(celltype_all, function(celltype){
  idx <- which(celltype_vec == celltype)
  colMeans(adt_mat[idx,,drop = F])
}))
adt_avg <- adt_avg[,1:10]

break_vec <- c(seq(min(adt_avg), 0, length.out = 11), seq(0, max(adt_avg), length.out = 11)[-1])
break_vec[1] <- break_vec[1]-1
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+1
base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
col_vec <- rev(grDevices::colorRampPalette(base_palette)(length(break_vec)-1))

png("../../../out/figures/main/bm_25antibody_adt-heatmap.png",
    height = 3500, width = 3500*ncol(adt_avg)/nrow(adt_avg), res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(adt_avg), 
      asp = (nrow(adt_avg)-1)/(ncol(adt_avg)-1),
      main = "", xlab = "", ylab = "",
      xaxt = "n", yaxt = "n", bty = "n",
      breaks = break_vec, col = col_vec)

# draw lines
n <- nrow(adt_avg)
halfspacing <- 1/(2*(n-1))
num_types <- sapply(celltype_list, length)
for(i in 1:(length(num_types)-1)){
  y_val <- 1-((2*halfspacing)*(sum(num_types[1:i])-1)+halfspacing)
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 5, col = "white")
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 4, lty = 2, col = "black")
}
graphics.off()

##############################

# consensus pca
Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "ADT"
mat_2 <- Matrix::t(bm[["ADT"]]@data[Seurat::VariableFeatures(object = bm),])

mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

set.seed(10)
consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = mat_1b, mat_2 = mat_2b,
                                           dims_1 = 1:30, dims_2 = 1:18,
                                           dims_consensus = 1:30,
                                           center_1 = T, center_2 = T,
                                           recenter_1 = F, recenter_2 = F,
                                           rescale_1 = F, rescale_2 = F,
                                           scale_1 = T, scale_2 = T,
                                           verbose = 1)
consensus_dimred <- consensus_pca$dimred_consensus

consensus_dimred <- .processing_func(consensus_dimred)
consensus_avg <- t(sapply(celltype_all, function(celltype){
  idx <- which(celltype_vec == celltype)
  colMeans(consensus_dimred[idx,,drop = F])
}))
consensus_avg <- consensus_avg[,1:10]

break_vec <- c(seq(min(consensus_avg), 0, length.out = 11), seq(0, max(consensus_avg), length.out = 11)[-1])
break_vec[1] <- break_vec[1]-1
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+1
base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
col_vec <- rev(grDevices::colorRampPalette(base_palette)(length(break_vec)-1))

png("../../../out/figures/main/bm_25antibody_consensusPCA-heatmap.png",
    height = 3500, width = 3500*ncol(consensus_avg)/nrow(consensus_avg), res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(consensus_avg), 
      asp = (nrow(consensus_avg)-1)/(ncol(consensus_avg)-1),
      main = "", xlab = "", ylab = "",
      xaxt = "n", yaxt = "n", bty = "n",
      breaks = break_vec, col = col_vec)

# draw lines
n <- nrow(consensus_avg)
halfspacing <- 1/(2*(n-1))
num_types <- sapply(celltype_list, length)
for(i in 1:(length(num_types)-1)){
  y_val <- 1-((2*halfspacing)*(sum(num_types[1:i])-1)+halfspacing)
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 5, col = "white")
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 4, lty = 2, col = "black")
}
graphics.off()

#############################

common_dimred <-  .processing_func(multiSVD_obj$tcca_obj$common_score)
common_avg <- t(sapply(celltype_all, function(celltype){
  idx <- which(celltype_vec == celltype)
  colMeans(common_dimred[idx,,drop = F])
}))
common_avg <- common_avg[,1:10]

break_vec <- c(seq(min(common_avg), 0, length.out = 11), seq(0, max(common_avg), length.out = 11)[-1])
break_vec[1] <- break_vec[1]-1
break_vec[length(break_vec)] <- break_vec[length(break_vec)]+1
base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
col_vec <- rev(grDevices::colorRampPalette(base_palette)(length(break_vec)-1))

png("../../../out/figures/main/bm_25antibody_tCCA-commmon-heatmap.png",
    height = 3500, width = 3500*ncol(common_avg)/nrow(common_avg), res = 500, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
image(tiltedCCA:::.rotate(common_avg), 
      asp = (nrow(common_avg)-1)/(ncol(common_avg)-1),
      main = "", xlab = "", ylab = "",
      xaxt = "n", yaxt = "n", bty = "n",
      breaks = break_vec, col = col_vec)

# draw lines
n <- nrow(common_avg)
halfspacing <- 1/(2*(n-1))
num_types <- sapply(celltype_list, length)
for(i in 1:(length(num_types)-1)){
  y_val <- 1-((2*halfspacing)*(sum(num_types[1:i])-1)+halfspacing)
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 5, col = "white")
  lines(x = c(0,1), y = rep(y_val, 2), lwd = 4, lty = 2, col = "black")
}
graphics.off()

