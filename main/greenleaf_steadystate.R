rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_tcca_RNA-ATAC.RData")
# load("../../../out/Writeup14p/Writeup14p_10x_greenleaf_tcca.RData")
source("greenleaf_colorPalette.R")

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
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(rna_common)
alignment_vec <- sapply(1:n, function(i){
  df <- data.frame(rna = rna_common[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec^val, rank(alignment_vec))
})
greenleaf$alignment <- alignment_vec^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(greenleaf$alignment), max(greenleaf$alignment), length.out = num_color)
color_vec <- sapply(greenleaf$alignment, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

save(color_vec, alignment_vec,
     scaling_grid, scaling_quality,
     date_of_run, session_info,
     file = "../../../out/main/10x_greenleaf_steadystate.RData")

png(paste0("../../../out/figures/main/10x_greenleaf_tcca_steadystate-full.png"),
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = c(4,4,4,0.5))
plot(x = greenleaf[["common_tcca"]]@cell.embeddings[,1],
     y = greenleaf[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     xlab = colnames(greenleaf[["common_tcca"]]@cell.embeddings)[1],
     ylab = colnames(greenleaf[["common_tcca"]]@cell.embeddings)[2],
     main = paste0("Human brain (10x, RNA+ATAC)\nAlignment between ATAC and common RNA"),
     xaxt = "n", yaxt = "n", bty = "n")
axis(side = 1)
axis(side = 2)
graphics.off()


png(paste0("../../../out/figures/main/10x_greenleaf_tcca_steadystate-full_cleaned.png"),
    height = 2000, width = 2000, units = "px", res = 500)
par(mar = c(0.5,0.5,0.5,0.5))
plot(x = greenleaf[["common_tcca"]]@cell.embeddings[,1],
     y = greenleaf[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = "",
     xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

###########

Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}
alignment_vec_alt <- sapply(1:n, function(i){
  df <- data.frame(rna = mat_1b[i,],
                   atac = atac_pred[i,])
  lm_res <- stats::lm(rna ~ atac, data = df)
  summary(lm_res)$r.squared
})

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec_alt^val, rank(alignment_vec_alt))
})
greenleaf$alignment_alt <- alignment_vec_alt^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(greenleaf$alignment_alt), max(greenleaf$alignment_alt), length.out = num_color)
color_vec <- sapply(greenleaf$alignment_alt, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

png(paste0("../../../out/figures/main/10x_greenleaf_tcca_steadystate_alt-full.png"),
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = c(4,4,4,0.5))
plot(x = greenleaf[["common_tcca"]]@cell.embeddings[,1],
     y = greenleaf[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     xlab = colnames(greenleaf[["common_tcca"]]@cell.embeddings)[1],
     ylab = colnames(greenleaf[["common_tcca"]]@cell.embeddings)[2],
     main = paste0("Human brain (10x, RNA+ATAC)\nAlignment between ATAC and RNA"),
     xaxt = "n", yaxt = "n", bty = "n")
axis(side = 1)
axis(side = 2)
graphics.off()

###########################3

greenleaf$alignment2 <- alignment_vec
keep_vec <- rep(1, ncol(greenleaf))
zz <- greenleaf[["common_tcca"]]@cell.embeddings
keep_vec[intersect(which(zz[,1] <= 0), which(zz[,2]>0))] <- 0
keep_vec[which(zz[,2] >= 10)] <- 0
keep_vec[which(zz[,2] <= -9)] <- 0
keep_vec[which(greenleaf$celltype %in% c("EC/Peric.", "IN1", "IN2", "IN3", "RG", "SP"))] <- 0
table(keep_vec)
greenleaf$keep <- keep_vec
greenleaf2 <- subset(greenleaf, keep == 1)
col_palette2 <- col_palette
names(col_palette2) <- as.character(1:length(col_palette2))
greenleaf2$tmp_label <- plyr::mapvalues(greenleaf2$celltype, from = names(col_palette), to = names(col_palette2))

plot1 <- Seurat::VlnPlot(greenleaf2, features = "alignment2", 
                         group.by = 'tmp_label', 
                         sort = "decreasing", pt.size = 0, cols = col_palette2) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_steadystate_violinplot_cleaned-subset.png"),
                plot1, device = "png", width = 4.5, height = 2, units = "in",
                dpi = 500)

######################################

target_genes <- c("GLI3", "EOMES", "NEUROG1", "OLIG2", "SOX3", "SOX10")
cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_mat <- rna_common[cell_idx[order_vec],]
atac_mat <- atac_pred[cell_idx[order_vec],]

rna_mat2 <- t(scale(t(rna_mat)))
atac_mat2 <- t(scale(t(atac_mat)))
diff_mat <- (rna_mat2 - atac_mat2)/rna_mat2
colnames(diff_mat) <- colnames(rna_common)
diff_mat <- abs(diff_mat[,target_genes])
diff_mat <- pmin(diff_mat, 15)
summary(diff_mat)

p <- ncol(diff_mat); n <- nrow(diff_mat)
pred_diff_mat <- sapply(1:p, function(j){
  print(j)
  
  tmp_df <- data.frame(y = diff_mat[,j], x = 1:n)
  reg_res <- npregfast::frfast(y ~ x, data = tmp_df)
  x_vec <- 1:n
  y_vec <- stats::predict(object = reg_res, newdata = data.frame(x = x_vec))
  y_vec$Estimation[,"Pred"]
})
colnames(pred_diff_mat) <- colnames(diff_mat)
summary(pred_diff_mat)

n <- nrow(pred_diff_mat)
for(i in 1:ncol(pred_diff_mat)){
  png(paste0("../../../out/figures/main/10x_greenleaf_steadystate_difference-",
             colnames(pred_diff_mat)[i], "_cleaned.png"),
      height = 500, width = 800, units = "px", res = 300)
  par(mar = c(2,2,0.5,0.5), bg = NA)
  plot(x = seq(0, 1, length.out = n),
       y = diff_mat[,i], pch = 16, col = rgb(0.75, 0.75, 0.75, 0.15),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
       ylim = c(0, 1.5*max(pred_diff_mat[,i])))
  lines(x = seq(0, 1, length.out = n),
        y = pred_diff_mat[,i], col = "white", lwd = 6)
  lines(x = seq(0, 1, length.out = n),
        y = pred_diff_mat[,i], col = "black", lwd = 4)
  axis(1, at = c(0, 0.5, 1), cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
  axis(2, label = F, cex.axis = 1.25, cex.lab = 1.25, lwd = 2)
  graphics.off()
}

