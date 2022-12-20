rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_tcca_RNA-geneActivity.RData")
load("../../../out/main/10x_greenleaf_developmentalGenes.RData")
# load("../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")
# load("../../out/main/10x_greenleaf_developmentalGenes.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

selected_variables <- selection_res$selected_variables

# for each gene, fit two np regression based on Lineage1 and Lineage2, and assign 
#  the pseudotime and branch for each gene
lineage_idx1 <- which(greenleaf$Lineage1 == 1)
lineage_idx2 <- which(greenleaf$Lineage2 == 1)
ordering1 <- order(greenleaf$pseudotime[lineage_idx1], decreasing = F)
ordering2 <- order(greenleaf$pseudotime[lineage_idx2], decreasing = F)
x_vec1 <- greenleaf$pseudotime[lineage_idx1[ordering1]]
x_vec2 <- greenleaf$pseudotime[lineage_idx2[ordering2]]
names(x_vec1) <- NULL
names(x_vec2) <- NULL

# construct the matrix that we'll be doing changepoints on
cell_vec <- colnames(greenleaf)[sort(unique(c(lineage_idx1, lineage_idx2)))]
tmp_mat <- greenleaf[["SCT"]]@data[sort(unique(c(greenleaf[["SCT"]]@var.features, selected_variables))),cell_vec]
tmp_mat <- t(as.matrix(tmp_mat))
gene_vec <- sort(unique(c(colnames(tmp_mat)[which(matrixStats::colSds(tmp_mat) >= 0.01)], selected_variables)))
tmp_mat <- tmp_mat[,gene_vec]
tmp_mat <- t(scale(tmp_mat))
svd_res <- irlba::irlba(tmp_mat, k = 50)
full_mat <- tcrossprod(eSVD2:::.mult_mat_vec(svd_res$u, svd_res$d), svd_res$v)
rownames(full_mat) <- gene_vec
colnames(full_mat) <- cell_vec

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

variable_summary <- sapply(1:length(selected_variables), function(i){
  print(paste0("Working on gene ", i, " out of ", length(selected_variables)))
  gene <- selected_variables[i]
  
  y_vec1 <- full_mat[gene,colnames(greenleaf)[lineage_idx1[ordering1]]]
  y_vec2 <- full_mat[gene,colnames(greenleaf)[lineage_idx2[ordering2]]]
  names(y_vec1) <- NULL
  names(y_vec2) <- NULL
  
  df1 <- data.frame(x = x_vec1,  y = y_vec1)
  np_res1 <- npregfast::frfast(y ~ x, data = df1)
  max_val1 <- max(np_res1$p[,1,1])
  pseudotime1 <- np_res1$x[which.max(np_res1$p[,1,1])]
  
  df2 <- data.frame(x = x_vec2,  y = y_vec2)
  np_res2 <- npregfast::frfast(y ~ x, data = df2)
  max_val2 <- max(np_res2$p[,1,1])
  pseudotime2 <- np_res2$x[which.max(np_res2$p[,1,1])]
  
  if(max_val1 > max_val2){
    return(c(Lineage = 1, pseudotime = pseudotime1, 
             max_val1 = max_val1, max_val2 = max_val2))
  } else {
    return(c(Lineage = 2, pseudotime = pseudotime2, 
             max_val1 = max_val1, max_val2 = max_val2))
  }
})
variable_summary <- t(variable_summary)
table(variable_summary[,"Lineage"])
rownames(variable_summary) <- selected_variables

sapply(unique(greenleaf$celltype), function(celltype){
  vec <- greenleaf$pseudotime[which(greenleaf$celltype == celltype)]
  vec <- vec[!is.na(vec)]
  stats::quantile(vec)
})

# hard set pseudotimes too early as lineage 1 for visualization purposes
branching_threshold <- 0.6
for(i in 1:nrow(variable_summary)){
  if(variable_summary[i,"pseudotime"] <= branching_threshold & variable_summary[i,"Lineage"] == 2){
    variable_summary[i,"Lineage"] <- 1
  }
}

# extract the relevant matrix
lineage_name_all <- colnames(greenleaf)[sort(unique(c(lineage_idx1, lineage_idx2)))]
heatmap_mat <- full_mat[selected_variables,lineage_name_all]
# heatmap_mat <- greenleaf[["SCT"]]@data[selected_variables,lineage_idx_all]
heatmap_mat <- t(as.matrix(heatmap_mat))

# order the genes in the matrix
tmp1 <- which(variable_summary[,"Lineage"] == 1)
gene_ordering1 <- rownames(variable_summary)[tmp1[order(variable_summary[tmp1,"pseudotime"], decreasing = F)]]
tmp2 <- which(variable_summary[,"Lineage"] == 2)
gene_ordering2 <- rownames(variable_summary)[tmp2[order(variable_summary[tmp2,"pseudotime"], decreasing = F)]]
gene_ordering <- c(gene_ordering1, gene_ordering2)
heatmap_mat <- heatmap_mat[,gene_ordering]

# split the matrix by cells per-lineage
cell_names <- rownames(heatmap_mat)
tmp <- which(colnames(greenleaf) %in% cell_names)
cell_names2 <- colnames(greenleaf)[tmp]
df_mat <- t(sapply(1:length(tmp), function(i){
  if(i %% floor(length(tmp)/10) == 0) cat('*')
  
  idx <- tmp[i]
  pseudotime <- greenleaf$pseudotime[idx]; names(pseudotime) <- NULL
  bool_lineage1 <- (greenleaf$Lineage1[idx] == 1)
  # if(pseudotime < branching_threshold | bool_lineage1){
  if(bool_lineage1){
    return(c(Lineage = 1, pseudotime = pseudotime))
  } else {
    return(c(Lineage = 2, pseudotime = pseudotime))
  }
}))
rownames(df_mat) <- cell_names2
tmp1 <- which(df_mat[,"Lineage"] == 1)
cell_ordering1 <- rownames(df_mat)[tmp1[order(df_mat[tmp1,"pseudotime"], decreasing = F)]]
tmp2 <- which(df_mat[,"Lineage"] == 2)
cell_ordering2 <- rownames(df_mat)[tmp2[order(df_mat[tmp2,"pseudotime"], decreasing = F)]]
cell_ordering <- c(cell_ordering1, cell_ordering2)

# split the matrices
scaling_val <- 1
heatmap_mat1 <- t(heatmap_mat[cell_ordering1,])
heatmap_mat1 <- sign(heatmap_mat1) * abs(heatmap_mat1)^scaling_val #^scaling_grid[which.max(scaling_quality)]
heatmap_mat2 <- t(heatmap_mat[cell_ordering2,])
heatmap_mat2 <- sign(heatmap_mat2) * abs(heatmap_mat2)^scaling_val # ^scaling_grid[which.max(scaling_quality)]

# # manual ordering
all(rownames(heatmap_mat1)[1:length(gene_ordering1)] == gene_ordering1)
all(rownames(heatmap_mat2)[(length(gene_ordering1)+1):nrow(heatmap_mat2)] == gene_ordering2)

heatmap_mat1b <- heatmap_mat1; heatmap_mat2b <- heatmap_mat2
gene_ordering1b <- gene_ordering1; gene_ordering2b <- gene_ordering2
heatmap_mat1b[(6:10),] <- heatmap_mat1b[c(10,6:9),]
heatmap_mat2b[(6:10),] <- heatmap_mat2b[c(10,6:9),]
rownames(heatmap_mat1b)[6:10] <- rownames(heatmap_mat1b)[c(10,6:9)]
rownames(heatmap_mat2b) <- rownames(heatmap_mat1b)
gene_ordering1b[6:10] <- gene_ordering1b[c(10,6:9)]

heatmap_mat1b[(11:14),] <- heatmap_mat1b[c(14,11:13),]
heatmap_mat2b[(11:14),] <- heatmap_mat2b[c(14,11:13),]
rownames(heatmap_mat1b)[11:14] <- rownames(heatmap_mat1b)[c(14,11:13)]
rownames(heatmap_mat2b) <- rownames(heatmap_mat1b)
gene_ordering1b[(11:14)] <- gene_ordering1b[c(14,11:13)]

gene_ordering2b <- c(gene_ordering2b[1], gene_ordering1b[c(38:40, 35)], gene_ordering2b[-1])
gene_ordering1b <- gene_ordering1b[c(1:34,36:37)]
heatmap_mat1b[(38:41),] <- heatmap_mat1b[c(41,38:40),]
heatmap_mat2b[(38:41),] <- heatmap_mat2b[c(41,38:40),]
rownames(heatmap_mat1b)[(38:41)] <- rownames(heatmap_mat1b)[c(41,38:40)]
rownames(heatmap_mat2b) <- rownames(heatmap_mat1b)
heatmap_mat1b[c(1:41),] <- heatmap_mat1b[c(1:34,36:41,35),]
heatmap_mat2b[c(1:41),] <- heatmap_mat2b[c(1:34,36:41,35),]
rownames(heatmap_mat1b)[c(1:41)] <- rownames(heatmap_mat1b)[c(1:34,36:41,35)]
rownames(heatmap_mat2b) <- rownames(heatmap_mat1b)
all(rownames(heatmap_mat1b)[1:length(gene_ordering1b)] == gene_ordering1b)
all(rownames(heatmap_mat2b)[(length(gene_ordering1b)+1):nrow(heatmap_mat2b)] == gene_ordering2b)

# set up the color palette
base_palette <- rev(RColorBrewer::brewer.pal(11, name = "RdYlBu"))
color_palette <- c(grDevices::colorRampPalette(base_palette[1:5])(5),
                   grDevices::colorRampPalette(base_palette[5:8])(70),
                   grDevices::colorRampPalette(base_palette[8:11])(25))
color_breaks <- seq(min(c(heatmap_mat1b, heatmap_mat2b)), 
                    max(c(heatmap_mat1b, heatmap_mat2b)), 
                    length.out = length(color_palette)+1)

height <- 1500
png(paste0("../../../out/figures/main/10x_greenleaf_developmentalGenes_heatmap1.png"),
    height = height, width = round(height*1.8/1.2*.8), units = "px", res = 500)
par(mar = c(0,0,0,0))
image(tiltedCCA:::.rotate(heatmap_mat1b), xaxs = "i", yaxs = "i",
      col = color_palette, breaks = color_breaks,
      main = "",
      xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")

num_down <- length(which(gene_ordering1b %in% rownames(heatmap_mat1b)))
num_total <- nrow(heatmap_mat1b)
y_val <- 1-num_down*(1/(num_total-1))+1/(2*(num_total-1))
lines(x = c(0,1), y = rep(y_val, 2), col = "white", lwd = 4)
lines(x = c(0,1), y = rep(y_val, 2), lty = 2, lwd = 3.5)
graphics.off()

png(paste0("../../../out/figures/main/10x_greenleaf_developmentalGenes_heatmap2.png"),
    height = height, width = round(height*1.8/1.2*.2), units = "px", res = 500)
par(mar = c(0,0,0,0))
image(tiltedCCA:::.rotate(heatmap_mat2b),
      col = color_palette, breaks = color_breaks,
      main = "",
      xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")

num_down <- length(which(gene_ordering1b %in% rownames(heatmap_mat2b)))
num_total <- nrow(heatmap_mat2b)
y_val <- 1-num_down*(1/(num_total-1))+1/(2*(num_total-1))
lines(x = c(0,1), y = rep(y_val, 2), col = "white", lwd = 4)
lines(x = c(0,1), y = rep(y_val, 2), lty = 2, lwd = 3.5)
graphics.off()

gene_ordering1c <- rownames(heatmap_mat1b)[which(rownames(heatmap_mat1b) %in% gene_ordering1b)]
gene_ordering2c <- rownames(heatmap_mat1b)[which(rownames(heatmap_mat1b) %in% gene_ordering2b)]
sink("../../../out/main/10x_greenleaf_developmentalGenes_heatmapOrder.txt")
for(i in 1:length(gene_ordering1c)){
  cat(paste0(i, ": ", gene_ordering1c[i]))
  cat("\n")
}
cat("\n")
cat("====")
cat("\n")
for(i in 1:length(gene_ordering2c)){
  cat(paste0(i+length(gene_ordering1c), ": ", gene_ordering2c[i]))
  cat("\n")
}
sink()

sink("../../../out/main/10x_greenleaf_developmentalGenes_heatmapOrder_gene_unformatted.txt")
for(x in rownames(heatmap_mat1b)){
  cat(paste0(x, "\n"))
}
sink()

sink("../../../out/main/10x_greenleaf_developmentalGenes_heatmapOrder_cell1_unformatted.txt")
for(x in colnames(heatmap_mat1b)){
  cat(paste0(x, "\n"))
}
sink()

sink("../../../out/main/10x_greenleaf_developmentalGenes_heatmapOrder_cell2_unformatted.txt")
for(x in colnames(heatmap_mat2b)){
  cat(paste0(x, "\n"))
}
sink()


##################################

# steady state plot
# load("../../../out/main/10x_greenleaf_steadystate.RData")
load("../../../out/main/10x_greenleaf_steadystate.RData")

alignment_vec1 <- alignment_vec[colnames(heatmap_mat1)]
df1 <- data.frame(y = alignment_vec1, x = 1:length(alignment_vec1))
np_res1 <- npregfast::frfast(y ~ x, data = df1)
col <- viridis::viridis(5)[3]

height <- 350
png(paste0("../../../out/figures/main/10x_greenleaf_steadystate_forHeatmap1.png"),
    height = height, width = height*2000/360, units = "px", res = 500)
par(mar = c(0.1,1,0.1,0), bg = NA)
y <- np_res1$p[,1,1]; n <- length(y)
ylim <- quantile(alignment_vec1, probs = c(0.1,0.95))
ylim[1] <- min(ylim[1], min(y)); ylim[2] <- max(ylim[2], max(y))
tmp <- alignment_vec1; tol <- .01
tmp[which(tmp <= ylim[1]+tol)] <- NA
tmp[which(tmp >= ylim[2]-tol)] <- NA
plot(x = seq(0, 1, length = length(alignment_vec1)),
     y = tmp, 
     ylim = ylim,
     col = rgb(0.5, 0.5, 0.5, 0.2), main = "", 
     xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "",
     pch = 16)
lines(x = seq(0, 1, length = n), y = y, lwd = 9, col = "white")
lines(x = seq(0, 1, length = n), y = y, lwd = 6, col = col)
axis(2, labels = F, lwd = 2)
graphics.off()

png(paste0("../../../out/figures/main/10x_greenleaf_steadystate_forHeatmap1-cleaned.png"),
    height = 900, width = 2500, units = "px", res = 500)
par(mar = c(2.5,2.5,0.5,0.5), bg = NA)
y <- np_res1$p[,1,1]; n <- length(y)
ylim <- quantile(alignment_vec1, probs = c(0.1,0.95))
ylim[1] <- min(ylim[1], min(y)); ylim[2] <- max(ylim[2], max(y))
tmp <- alignment_vec1; tol <- .01
tmp[which(tmp <= ylim[1]+tol)] <- NA
tmp[which(tmp >= ylim[2]-tol)] <- NA
plot(x = seq(0, 1, length = length(alignment_vec1)),
     y = tmp, 
     ylim = ylim,
     col = rgb(0.5, 0.5, 0.5, 0.2), main = "", 
     bty = "n", xlab = "", ylab = "",
     pch = 16)
lines(x = seq(0, 1, length = n), y = y, lwd = 9, col = "white")
lines(x = seq(0, 1, length = n), y = y, lwd = 6, col = col)
graphics.off()

# the second lineage is tricker since 1) we want to fit the curve on all the points in
#  the second lineage, and 2) we only plot the predictions for the relevant cells in heatmap_mat2,
#  but np_res2 outputs the regression values on a grid
lineage_idx2 <- which(greenleaf$Lineage2 == 1)
lineage_idx_reordered <- lineage_idx2[order(greenleaf$pseudotime[lineage_idx2])]
alignment_vec2 <- alignment_vec[lineage_idx_reordered]
names(alignment_vec2) <- colnames(greenleaf)[lineage_idx_reordered]
x_full <- 1:length(alignment_vec2)
names(x_full) <- colnames(greenleaf)[lineage_idx_reordered]
df2 <- data.frame(y = alignment_vec2, x = x_full)
np_res2 <- npregfast::frfast(y ~ x, data = df2)

quantile(np_res2$p[,1,1])
alignment_vec2_subset <- alignment_vec2[colnames(heatmap_mat2)]
x_subset <- x_full[colnames(heatmap_mat2)]
predicted_y <- stats::predict(np_res2, newdata = data.frame(x = x_subset))
predicted_y <- predicted_y$Estimation[,"Pred"]

png(paste0("../../../out/figures/main/10x_greenleaf_steadystate_forHeatmap2.png"),
    height = height, width = height*1900/360/.8*.2, units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0), bg = NA)
n <- length(x_subset)
tmp <- alignment_vec2_subset; tol <- .01
tmp[which(tmp <= ylim[1]+tol)] <- NA
tmp[which(tmp >= ylim[2]-tol)] <- NA
plot(x = seq(0, 1, length = n),
     y = tmp, 
     col = rgb(0.5, 0.5, 0.5, 0.2), main = "", 
     ylim = ylim,
     xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "",
     pch = 16)
lines(x = seq(0, 1, length = n), y = predicted_y, lwd = 9, col = "white")
lines(x = seq(0, 1, length = n), y = predicted_y, lwd = 6, col = col)
graphics.off()

#####################################################

source("../tiltedCCA_analysis/main/greenleaf_colorPalette.R")
celltype_vec1 <- sapply(colnames(heatmap_mat1), function(cell_name){
  greenleaf$celltype[cell_name]
})
celltype_vec1 <- as.factor(celltype_vec1)
uniq_names <- levels(celltype_vec1)
celltype_vec1 <- as.numeric(celltype_vec1)
color_vec <- col_palette[uniq_names]

png(paste0("../../out/figures/main/10x_greenleaf_celltypes_forHeatmap1.png"),
    height = 350, width = height*1700/150, units = "px", res = 500)
par(mar = c(0,0,0,0))
image(as.matrix(celltype_vec1, ncol = length(celltype_vec1), nrow = 1),
      breaks = seq(.5, length(uniq_names)+.5, by = 1),
      col = color_vec)
graphics.off()


celltype_vec2 <- sapply(colnames(heatmap_mat2), function(cell_name){
  greenleaf$celltype[cell_name]
})
celltype_vec2 <- as.factor(celltype_vec2)
uniq_names <- levels(celltype_vec2)
celltype_vec2 <- as.numeric(celltype_vec2)
color_vec <- col_palette[uniq_names]
png(paste0("../../out/figures/main/10x_greenleaf_celltypes_forHeatmap2.png"),
    height = 350, width = height*1700/150/.8*.2, units = "px", res = 500)
par(mar = c(0,0,0,0))
image(as.matrix(celltype_vec2, ncol = length(celltype_vec2), nrow = 1),
      breaks = seq(.5, length(uniq_names)+.5, by = 1),
      col = color_vec)
graphics.off()

