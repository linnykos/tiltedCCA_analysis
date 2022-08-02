rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_reik_tiltedcca.RData")
source("reik_colorPalette.R")

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
round(quantile(alignment_vec), 2)

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec^val, rank(alignment_vec))
})
reik$alignment <- alignment_vec^(scaling_grid[which.max(scaling_quality)])
reik$alignment2 <- alignment_vec
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(reik$alignment), max(reik$alignment), length.out = num_color)
color_vec <- sapply(reik$alignment, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

save(color_vec, alignment_vec,
     scaling_grid, scaling_quality,
     date_of_run, session_info,
     file = "../../../out/main/10x_reik_steadystate.RData")

png(paste0("../../../out/figures/main/10x_reik_tcca_steadystate-full.png"),
    height = 3000, width = 3000, units = "px", res = 500)
par(mar = c(4,4,4,0.5))
plot(x = reik[["common_tcca"]]@cell.embeddings[,1],
     y = reik[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     xlab = colnames(reik[["common_tcca"]]@cell.embeddings)[1],
     ylab = colnames(reik[["common_tcca"]]@cell.embeddings)[2],
     main = paste0("Mouse embryo (10x, RNA+ATAC)\nAlignment between ATAC and common RNA"),
     xaxt = "n", yaxt = "n", bty = "n")
axis(side = 1)
axis(side = 2)
graphics.off()


plot1 <- Seurat::VlnPlot(reik, features = "alignment2", 
                         group.by = 'celltype', 
                         sort = TRUE, pt.size = 0.1, cols = col_palette) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_steadystate_violinplot.png"),
                plot1, device = "png", width = 20, height = 5, units = "in")

############################

idx <- intersect(which(!is.na(reik$celltype)), which(!reik$celltype == "ExE_ectoderm"))

scaling_grid <- seq(0.1, 10, length.out = 100)
scaling_quality <- sapply(scaling_grid, function(val){
  stats::cor(alignment_vec[idx]^val, rank(alignment_vec[idx]))
})

custom_alignment_vec <- alignment_vec[idx]^(scaling_grid[which.max(scaling_quality)])
num_color <- 100
color_palette <- viridis::viridis(num_color)
color_breaks <- seq(min(custom_alignment_vec), max(custom_alignment_vec), length.out = num_color)
color_vec <- sapply(custom_alignment_vec, function(val){
  color_palette[which.min(abs(color_breaks - val))]
})

png(paste0("../../../out/figures/main/10x_reik_tcca_steadystate-full_cleaned.png"),
    height = 2000, width = 2000, units = "px", res = 500)
par(mar = c(0.5,0.5,0.5,0.5))
plot(x = reik[["common_tcca"]]@cell.embeddings[idx,1],
     y = reik[["common_tcca"]]@cell.embeddings[idx,2],
     col = color_vec, pch = 16,
     main = "",
     xaxt = "n", yaxt = "n", bty = "n")
graphics.off()

## 

idx <- intersect(which(!is.na(reik$celltype)), which(!reik$celltype == "ExE_ectoderm"))
keep_vec <- rep(0, ncol(reik))
keep_vec[idx] <- 1
reik$keep <- keep_vec
reik2 <- subset(reik, keep == 1)

custom_celltype <- rep("NA", ncol(reik2))
celltypes_keep <- c("Cardiomyocytes", "ExE_endoderm", "Parietal_endoderm",
                    "Erythroid3", "Blood_progenitors_1", "Primitive_Streak",
                    "Anterior_Primitive_Streak", "Nascent_mesoderm",
                    "Endothelium")
color_vec <- rep(rgb(0.8,0.8,0.8), ncol(reik2))
for(celltype in celltypes_keep){
  color_vec[which(reik2$celltype == celltype)] <- col_palette[celltype]
}

png(paste0("../../../out/figures/main/10x_reik_tcca-umap_common_cleaned-smaller.png"),
    height = 2000, width = 2000, units = "px", res = 500)
par(mar = c(0.5,0.5,0.5,0.5))
plot(x = reik2[["common_tcca"]]@cell.embeddings[,1],
     y = reik2[["common_tcca"]]@cell.embeddings[,2],
     col = color_vec, pch = 16,
     main = "",
     xaxt = "n", yaxt = "n", bty = "n")
idx <- which(reik2$celltype %in% celltypes_keep)
# points(x = reik2[["common_tcca"]]@cell.embeddings[idx,1],
#        y = reik2[["common_tcca"]]@cell.embeddings[idx,2],
#        col = "white", pch = 16, cex = 2)
points(x = reik2[["common_tcca"]]@cell.embeddings[idx,1],
       y = reik2[["common_tcca"]]@cell.embeddings[idx,2],
       col = color_vec[idx], pch = 16, cex = 1)
graphics.off()

#####################################

reik$alignment2 <- alignment_vec
keep_vec <- rep(1, ncol(reik))
keep_vec[is.na(reik$celltype)] <- 0
celltypes_keep <- c("Cardiomyocytes", "ExE_endoderm", "Parietal_endoderm",
                    "Erythroid3", "Blood_progenitors_1", "Primitive_Streak",
                    "Anterior_Primitive_Streak", "Nascent_mesoderm",
                    "Endothelium")
keep_vec[!reik$celltype %in% celltypes_keep] <- 0
reik$keep <- keep_vec
reik2 <- subset(reik, keep == 1)
col_palette2 <- col_palette
names(col_palette2) <- as.character(1:length(col_palette2))
reik2$tmp_label <- plyr::mapvalues(reik2$celltype, from = names(col_palette), to = names(col_palette2))

plot1 <- Seurat::VlnPlot(reik2, features = "alignment2", 
                         group.by = 'tmp_label', 
                         sort = "decreasing", pt.size = 0, cols = col_palette2) + Seurat::NoLegend()
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_reik_steadystate_violinplot_cleaned-subset.png"),
                plot1, device = "png", width = 5.6, height = 2, units = "in",
                dpi = 500)
