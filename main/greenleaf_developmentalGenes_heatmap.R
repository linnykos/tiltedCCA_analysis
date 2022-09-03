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

variable_summary <- sapply(1:length(selected_variables), function(i){
  print(paste0("Working on gene ", i, " out of ", length(selected_variables)))
  gene <- selected_variables[i]
  
  y_vec1 <- greenleaf[["SCT"]]@data[gene,lineage_idx1[ordering1]]
  y_vec2 <- greenleaf[["SCT"]]@data[gene,lineage_idx2[ordering2]]
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
    return(c(Lineage = 1, pseudotime = pseudotime1))
  } else {
    return(c(Lineage = 2, pseudotime = pseudotime2))
  }
})
variable_summary <- t(variable_summary)
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
lineage_idx_all <- sort(unique(c(lineage_idx1, lineage_idx2)))
heatmap_mat <- greenleaf[["SCT"]]@data[selected_variables,lineage_idx_all]
heatmap_mat <- as.matrix(heatmap_mat)

# normalize the cells in the matrix
heatmap_mat <- scale(t(heatmap_mat))

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

# threshold the values
heatmap_mat_threshold <- heatmap_mat
for(j in 1:ncol(heatmap_mat_threshold)){
  max_val <- quantile(heatmap_mat_threshold[,j], probs = .9)
  min_val <- quantile(heatmap_mat_threshold[,j], probs = .1)
  if(abs(max_val - min_val) <= .1){
    max_val <- quantile(heatmap_mat_threshold[,j], probs = .95)
    min_val <- quantile(heatmap_mat_threshold[,j], probs = .05)
  }
  heatmap_mat_threshold[,j] <- pmax(pmin(heatmap_mat_threshold[,j], max_val), min_val)
}

heatmap_mat_threshold <- heatmap_mat_threshold[,which(apply(heatmap_mat_threshold, 2, sd) > 0.1)]

# split the matrices
scaling_val <- 0.5
heatmap_mat1 <- t(heatmap_mat_threshold[cell_ordering1,])
heatmap_mat1 <- sign(heatmap_mat1) * abs(heatmap_mat1)^scaling_val #^scaling_grid[which.max(scaling_quality)]
heatmap_mat2 <- t(heatmap_mat_threshold[cell_ordering2,])
heatmap_mat2 <- sign(heatmap_mat2) * abs(heatmap_mat2)^scaling_val # ^scaling_grid[which.max(scaling_quality)]

# set up the color palette
num_color <- 100
base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_palette <- rev(grDevices::colorRampPalette(base_palette)(100))
color_breaks <- seq(min(c(heatmap_mat1, heatmap_mat2)), 
                    max(c(heatmap_mat1, heatmap_mat2)), 
                    length.out = num_color+1)

height <- 1500
png(paste0("../../../out/figures/main/10x_greenleaf_developmentalGenes_heatmap1.png"),
    height = height, width = round(height*1.8/1.2*.8), units = "px", res = 500)
par(mar = c(0,0,0,0))
image(tiltedCCA:::.rotate(heatmap_mat1), xaxs = "i", yaxs = "i",
      col = color_palette, breaks = color_breaks,
      main = "",
      xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")

num_down <- length(which(gene_ordering1 %in% rownames(heatmap_mat1)))
num_total <- nrow(heatmap_mat1)
y_val <- 1-num_down*(1/(num_total-1))+1/(2*(num_total-1))
lines(x = c(0,1), y = rep(y_val, 2), col = "white", lwd = 4)
lines(x = c(0,1), y = rep(y_val, 2), lty = 2, lwd = 3.5)
graphics.off()

png(paste0("../../../out/figures/main/10x_greenleaf_developmentalGenes_heatmap2.png"),
    height = height, width = round(height*1.8/1.2*.2), units = "px", res = 500)
par(mar = c(0,0,0,0))
image(tiltedCCA:::.rotate(heatmap_mat2),
      col = color_palette, breaks = color_breaks,
      main = "",
      xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "")

num_down <- length(which(gene_ordering1 %in% rownames(heatmap_mat1)))
num_total <- nrow(heatmap_mat1)
y_val <- 1-num_down*(1/(num_total-1))+1/(2*(num_total-1))
lines(x = c(0,1), y = rep(y_val, 2), col = "white", lwd = 4)
lines(x = c(0,1), y = rep(y_val, 2), lty = 2, lwd = 3.5)
graphics.off()

gene_ordering1b <- gene_ordering1[which(gene_ordering1 %in% rownames(heatmap_mat1))]
gene_ordering2b <- gene_ordering2[which(gene_ordering2 %in% rownames(heatmap_mat1))]
sink("../../../out/main/10x_greenleaf_developmentalGenes_heatmapOrder.txt")
for(i in 1:length(gene_ordering1b)){
  cat(paste0(i, ": ", gene_ordering1b[i]))
  cat("\n")
}
cat("\n")
cat("====")
cat("\n")
for(i in 1:length(gene_ordering2b)){
  cat(paste0(i+length(gene_ordering1b), ": ", gene_ordering2b[i]))
  cat("\n")
}
sink()

##################################

# steady state plot
# load("../../../out/main/10x_greenleaf_steadystate.RData")
load("../../out/main/10x_greenleaf_steadystate.RData")

names(alignment_vec) <- colnames(greenleaf)
alignment_vec1 <- alignment_vec[colnames(heatmap_mat1)]
df1 <- data.frame(y = alignment_vec1, x = 1:length(alignment_vec1))
np_res1 <- npregfast::frfast(y ~ x, data = df1)
col <- viridis::viridis(5)[3]

height <- 350
png(paste0("../../out/figures/main/10x_greenleaf_steadystate_forHeatmap1.png"),
    height = height, width = height*1800/360, units = "px", res = 500)
par(mar = c(0.1,1,0.1,0))
y <- np_res1$p[,1,1]; n <- length(y)
ylim <- quantile(alignment_vec1, probs = c(0.1,0.95))
ylim[1] <- min(ylim[1], min(y)); ylim[2] <- max(ylim[2], max(y))
tmp <- alignment_vec1; tol <- .01
tmp[which(tmp <= ylim[1]+tol)] <- NA
tmp[which(tmp >= ylim[2]-tol)] <- NA
plot(x = seq(0, 1, length = length(alignment_vec1)),
     y = tmp, 
     ylim = ylim,
     col = rgb(0.5, 0.5, 0.5, 0.1), main = "", 
     xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "",
     pch = 16)
lines(x = seq(0, 1, length = n), y = y, lwd = 9, col = "white")
lines(x = seq(0, 1, length = n), y = y, lwd = 6, col = col)
axis(2, labels = F, lwd = 2)
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

png(paste0("../../out/figures/main/10x_greenleaf_steadystate_forHeatmap2.png"),
    height = height, width = height*1700/360/.8*.2, units = "px", res = 500)
par(mar = c(0.1,0.1,0.1,0))
n <- length(x_subset)
tmp <- alignment_vec2_subset; tol <- .01
tmp[which(tmp <= ylim[1]+tol)] <- NA
tmp[which(tmp >= ylim[2]-tol)] <- NA
plot(x = seq(0, 1, length = n),
     y = tmp, 
     col = rgb(0.5, 0.5, 0.5, 0.1), main = "", 
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

