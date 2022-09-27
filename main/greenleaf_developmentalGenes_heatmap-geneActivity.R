rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_tcca_RNA-geneActivity.RData")
load("../../../out/main/10x_greenleaf_developmentalGenes.RData")

gene_ordering <- read.csv("../../../out/main/10x_greenleaf_developmentalGenes_heatmapOrder_gene_unformatted.txt",
                          header = F)
cell1_ordering <- read.csv("../../../out/main/10x_greenleaf_developmentalGenes_heatmapOrder_cell1_unformatted.txt",
                          header = F)
cell2_ordering <- read.csv("../../../out/main/10x_greenleaf_developmentalGenes_heatmapOrder_cell2_unformatted.txt",
                          header = F)

#######

# construct the matrix that we'll be doing changepoints on
cell_vec <- unique(c(cell1_ordering[,1], cell2_ordering[,1]))
# tmp_mat <- greenleaf[["customGAct"]]@data[greenleaf[["customGAct"]]@var.features,cell_vec]
tmp_mat <- t(multiSVD_obj$common_mat_2)
full_mat <- tmp_mat[greenleaf[["customGAct"]]@var.features,cell_vec]
# tmp_mat <- t(as.matrix(tmp_mat))
# gene_vec <- sort(unique(c(colnames(tmp_mat)[which(matrixStats::colSds(tmp_mat) >= 0.01)], paste0("ATAC-", gene_ordering[,1]))))
# tmp_mat <- tmp_mat[,gene_vec]
# tmp_mat <- t(scale(tmp_mat))
# svd_res <- irlba::irlba(tmp_mat, k = 50)
# full_mat <- tcrossprod(eSVD2:::.mult_mat_vec(svd_res$u, svd_res$d), svd_res$v)
# rownames(full_mat) <- gene_vec
# colnames(full_mat) <- cell_vec

# threshold the values
for(j in 1:nrow(full_mat)){
  max_val <- quantile(full_mat[j,], probs = .9)
  min_val <- quantile(full_mat[j,], probs = .1)
  if(abs(max_val - min_val) <= .1){
    max_val <- quantile(full_mat[j,], probs = .95)
    min_val <- quantile(full_mat[j,], probs = .05)
  }
  full_mat[j,] <- pmax(pmin(full_mat[j,], max_val), min_val)
  full_mat[j,] <- (full_mat[j,] - min_val)/(max_val - min_val)
}

heatmap_mat1b <- full_mat[paste0("ATAC-", gene_ordering[,1]), cell1_ordering[,1]]
heatmap_mat2b <- full_mat[paste0("ATAC-", gene_ordering[,1]), cell2_ordering[,1]]

# set up the color palette
base_palette <- rev(RColorBrewer::brewer.pal(11, name = "RdYlBu"))
color_palette <- c(grDevices::colorRampPalette(base_palette[1:5])(20),
                   grDevices::colorRampPalette(base_palette[5:8])(70),
                   grDevices::colorRampPalette(base_palette[8:11])(10))
color_breaks <- seq(min(c(heatmap_mat1b, heatmap_mat2b)), 
                    max(c(heatmap_mat1b, heatmap_mat2b)), 
                    length.out = length(color_palette)+1)

height <- 1500
png(paste0("../../../out/figures/main/10x_greenleaf_developmentalGenes-geneActivity_heatmap1.png"),
    height = height, width = round(height*1.8/1.2*.8), units = "px", res = 500)
par(mar = c(0,0,0,0))
image(tiltedCCA:::.rotate(heatmap_mat1b), xaxs = "i", yaxs = "i",
      col = color_palette, breaks = color_breaks,
      main = "",
      xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
graphics.off()

png(paste0("../../../out/figures/main/10x_greenleaf_developmentalGenes-geneActivity_heatmap2.png"),
    height = height, width = round(height*1.8/1.2*.2), units = "px", res = 500)
par(mar = c(0,0,0,0))
image(tiltedCCA:::.rotate(heatmap_mat2b),
      col = color_palette, breaks = color_breaks,
      main = "",
      xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")
graphics.off()