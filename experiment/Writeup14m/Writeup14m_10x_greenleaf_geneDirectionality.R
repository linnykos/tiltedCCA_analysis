rm(list=ls())
load("../../../../out/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_tcca.RData")
library(Seurat)
source("../Writeup14l/slingshot.R")

set.seed(10)
greenleaf <- Seurat::FindMultiModalNeighbors(greenleaf, reduction.list = list("pca", "lsi"), 
                                             dims.list = list(1:50, 2:50))
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf,
                                  graph.name = "wsnn", algorithm = 3, 
                                  resolution = 2)
lineage_order <- as.character(c(11,15,8,3,14,0,1,2,6,12))

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nSeurat WNN clusters, on Common embedding"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_activity_seuratWNNcluster.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

cell_names <- colnames(greenleaf)[which(greenleaf$seurat_clusters %in% lineage_order)]
plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         cells.highlight = cell_names)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nSeurat WNN clusters (subset)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_activity_seuratWNNcluster_subset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")


#########################

multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
dimred_1 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, 
                                     apply_postDimred = T, what = "common_mat")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred_2 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, what = "common_mat")
dimred <- cbind(dimred_1, dimred_2)
svd_tmp <- irlba::irlba(dimred, nv = 25)
dimred <- tiltedCCA:::.mult_mat_vec(svd_tmp$u[,1:10], svd_tmp$d[1:10])

cluster_vec <- greenleaf$seurat_clusters[which(greenleaf$seurat_clusters %in% lineage_order)]
dimred <- dimred[which(greenleaf$seurat_clusters %in% lineage_order),]

initial_fit <- .initial_curve_fit(cluster_vec = cluster_vec,
                                  dimred = dimred,
                                  lineage_order = lineage_order)
pseudotime_vec <- .extract_pseudotime(dimred = dimred,
                                      initial_fit = initial_fit,
                                      stretch = 2)
n <- ncol(greenleaf)
pseudotime_full <- rep(NA, n)
pseudotime_full[which(greenleaf$seurat_clusters %in% lineage_order)] <- pseudotime_vec
pseudotime_full2 <- rep(NA, n)
pseudotime_full2[which(greenleaf$seurat_clusters %in% lineage_order)] <- rank(pseudotime_vec)

greenleaf$pseudotime <- pseudotime_full
greenleaf$pseudotime_rank <- pseudotime_full2
sapply(lineage_order, function(clust){
  idx <- which(greenleaf$seurat_clusters == clust)
  if(length(idx) == 0) return(NA)
  mean(greenleaf$pseudotime_rank[idx], na.rm = T)
})
plot1 <- Seurat::FeaturePlot(greenleaf, feature = "pseudotime",
                             reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nSlingshot pseudotime via Common embedding"))
plot2 <- Seurat::FeaturePlot(greenleaf, feature = "pseudotime_rank",
                             reduction = "common_tcca")
plot2 <- plot2 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nPseudotime (Ranked)"))

plot3 <- cowplot::plot_grid(plot1, plot2, ncol = 2)
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_activity_pseudotime.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")

################################

multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
dimred_1 <- tiltedCCA:::.get_postDimred(multiSVD_obj, averaging_mat = NULL)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred_2 <- tiltedCCA:::.get_postDimred(multiSVD_obj, averaging_mat = NULL)
dimred <- cbind(dimred_1, dimred_2)
param <- multiSVD_obj$param

snn_mat <- tiltedCCA:::.form_snn_mat(mat = dimred, 
                                     num_neigh = param$snn_num_neigh,
                                     bool_cosine = param$snn_bool_cosine,
                                     bool_intersect = param$snn_bool_intersect,
                                     min_deg = param$snn_min_deg)

# rna_common <- multiSVD_obj$common_mat_1
# rna_distinct <- multiSVD_obj$distinct_mat_1
# n <- nrow(rna_common)
# alignment_vec <- sapply(1:n, function(i){
#   df <- data.frame(common = rna_common[i,],
#                    distinct = rna_distinct[i,])
#   lm_res <- stats::lm(distinct ~ common, data = df)
#   1-tiltedCCA:::.l2norm(lm_res$residuals)/tiltedCCA:::.l2norm(rna_distinct[i,])
# })
# quantile(alignment_vec)
# 
# alignment_vec_smoothed <- sapply(1:n, function(i){
#   idx <- c(tiltedCCA:::.nonzero_col(snn_mat, col_idx = i, bool_value = F), i)
#   mean(alignment_vec[idx])
# })
# 
# greenleaf$alignment <- alignment_vec_smoothed

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

greenleaf$alignment <- alignment_vec
plot1 <-Seurat::FeaturePlot(greenleaf, feature = "alignment",
                            reduction = "common_tcca")
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nDistinct resid. (predict using Common, Smoothed)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14m/tmp.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

################################

rna_common <- multiSVD_obj$common_mat_1   
atac_everything <- multiSVD_obj$common_mat_2 + multiSVD_obj$distinct_mat_2

gene_names <- sapply(colnames(atac_everything), function(x){
  tmp <- strsplit(x, split="-")[[1]]
  paste0(tmp[2:length(tmp)], collapse="-")
})
names(gene_names) <- NULL
all(gene_names %in% colnames(rna_common))
rna_common <- rna_common[,gene_names]
all(colnames(rna_common) == colnames(gene_names))

cell_idx <- which(greenleaf$seurat_clusters %in% lineage_order)
rna_common <- rna_common[cell_idx,]
atac_everything <- atac_everything[cell_idx,]
rna_common <- scale(rna_common)
atac_everything <- scale(atac_everything)

p <- ncol(rna_common)
cor_vec <- sapply(1:p, function(j){
  vec1 <- rna_common[,j]
  vec2 <- atac_everything[,j]
  stats::cor(vec1, vec2)
})
quantile(cor_vec)

idx <- order(cor_vec, decreasing = T)[1:6]
col_pal <- grDevices::colorRampPalette(c('lightgrey', 'blue'))(100)
seq_vec <- seq(1, length(cell_idx), length.out = 100)
tmp <- greenleaf$pseudotime[cell_idx]
tmp <- tmp + runif(length(tmp), min = 0, max = 0.001)
tmp <- rank(tmp)
col_vec <- sapply(tmp, function(val){
  col_pal[which.min(abs(val - seq_vec))]
})
shuff_idx <- sample(1:nrow(rna_common))
rna_common <- rna_common[shuff_idx,]; atac_everything <- atac_everything[shuff_idx,]
col_vec <- col_vec[shuff_idx]
png("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_activity_everythingATAC2commonRNA_individualgenes.png",
    height = 2000, width = 3000, units = "px", res = 300)
par(mfrow = c(2,3))
for(j in idx){
  plot(atac_everything[,j], rna_common[,j],
       xlab = paste0("Everything ATAC, Gene activity for ", colnames(rna_common)[j]),
       ylab = paste0("Common RNA for ", colnames(rna_common)[j]),
       col = col_vec, pch = 16, cex = 0.5
  )
}
graphics.off()

########################

set.seed(10)
celltype_vec <- greenleaf$celltype[cell_idx]
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
pseudotime_vec <- pseudotime_vec + runif(length(pseudotime_vec), min = 0, max = 0.001)
order_vec <- order(pseudotime_vec, decreasing = F)
alignment_vec <- greenleaf$alignment[cell_idx]

rna_common2 <- rna_common[order_vec,]
atac_everything2 <- atac_everything[order_vec,]
celltype_vec <- celltype_vec[order_vec]
alignment_vec <- alignment_vec[order_vec]

p <- ncol(rna_common2); n <- nrow(rna_common2)
rna_unimodal_fit_mat <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  UniIsoRegression::reg_1d(rna_common2[,j],
                           w_vec = rep(1, n),
                           metric = 2,
                           unimodal = T)
})
atac_unimodal_fit_mat <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  UniIsoRegression::reg_1d(atac_everything2[,j],
                           w_vec = rep(1, n),
                           metric = 2,
                           unimodal = T)
})
rna_quality <- sapply(1:p, function(j){
  unimodal_mse <- mean((rna_unimodal_fit_mat[,j] - rna_common2[,j])^2)
  flat_mse <- sd(rna_common2[,j])
  abs(unimodal_mse - flat_mse)/unimodal_mse
})
atac_quality <- sapply(1:p, function(j){
  unimodal_mse <- mean((atac_unimodal_fit_mat[,j] - atac_everything2[,j])^2)
  flat_mse <- sd(atac_everything2[,j])
  abs(unimodal_mse - flat_mse)/unimodal_mse
})

num_considered <- 100
# idx1 <- order(cor_vec, decreasing = T)[1:num_considered]
idx2 <- order(rna_quality, decreasing = T)[1:num_considered]
idx3 <- order(atac_quality, decreasing = T)[1:num_considered]
# idx <- intersect(idx1, intersect(idx2, idx3))
idx <- intersect(idx2, idx3)
length(idx)
col_vec <- grDevices::colorRampPalette(c('lightgrey', 'blue'))(nrow(rna_common2))
png("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_activity_everythingATAC2commonRNA_individualgenes_unimodal_upward.png",
    height = 2000, width = 3000, units = "px", res = 300)
par(mfrow = c(3,5), mar = c(4,4,4,0.5))
for(j in idx[1:15]){
  plot(atac_everything2[,j], rna_common2[,j],
       main = paste0("Correlation: ", round(cor_vec[j], 2)),
       xlab = paste0("Gene activity for ", colnames(rna_common)[j]),
       ylab = paste0("Common RNA for ", colnames(rna_common2)[j]),
       col = col_vec, pch = 16, cex = 0.5
  )
}
graphics.off()
png("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_activity_everythingATAC2commonRNA_individualgenes_unimodal_upward2.png",
    height = 2000, width = 3000, units = "px", res = 300)
par(mfrow = c(3,5), mar = c(4,4,4,0.5))
for(j in idx[1:15]){
  plot(NA, 
       xlim = c(1,n),
       ylim = range(c(atac_everything2[,j], rna_common2[,j])),
       xlab = "Rank (via Slingshot)",
       ylab = paste0("Expression for ", colnames(rna_common2)[j])
  )
  
  x_vec <- c(1:n, 1:n)
  y_vec <- c(rna_common2[,j], atac_everything2[,j])
  col_vec <- c(rep(rgb(0.5, 0.5, 0.5, 0.5), n),
               rep(rgb(0.75, 0, 0, 0.5), n))
  shuff_idx <- sample(1:length(x_vec))
  points(x = x_vec[shuff_idx], y = y_vec[shuff_idx], 
         col = col_vec[shuff_idx], pch = 16, cex = 0.5)
  lines(x = 1:n, y = rna_unimodal_fit_mat[,j],
        col = "white", lwd = 5)
  lines(x = 1:n, y = rna_unimodal_fit_mat[,j],
        col = "black", lwd = 3)
  lines(x = 1:n, y = atac_unimodal_fit_mat[,j],
        col = "white", lwd = 5)
  lines(x = 1:n, y = atac_unimodal_fit_mat[,j],
        col = 2, lwd = 3)
}
graphics.off()

col_pal <- grDevices::colorRampPalette(c('lightgrey', 'blue'))(100)
seq_vec <- seq(min(alignment_vec), max(alignment_vec), length.out = 100)
col_vec <- sapply(alignment_vec, function(val){
  col_pal[which.min(abs(val - seq_vec))]
})
png("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_activity_everythingATAC2commonRNA_individualgenes_unimodal_upward_alignment.png",
    height = 2000, width = 3000, units = "px", res = 300)
par(mfrow = c(3,5), mar = c(4,4,4,0.5))
for(j in idx[1:15]){
  plot(atac_everything2[,j], rna_common2[,j],
       xlab = paste0("Gene activity for ", colnames(rna_common)[j]),
       ylab = paste0("Common RNA for ", colnames(rna_common2)[j]),
       col = col_vec, pch = 16, cex = 0.5
  )
}
graphics.off()

####################################
####################################
####################################

p <- ncol(rna_common2); n <- nrow(rna_common2)
rna_unimodal_fit_mat <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  UniIsoRegression::reg_1d(rna_common2[,j],
                           w_vec = rep(1, n),
                           metric = 2,
                           unimodal = F,
                           decreasing = T)
})
atac_unimodal_fit_mat <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  UniIsoRegression::reg_1d(atac_everything2[,j],
                           w_vec = rep(1, n),
                           metric = 2,
                           unimodal = F,
                           decreasing = T)
})
rna_quality <- sapply(1:p, function(j){
  unimodal_mse <- mean((rna_unimodal_fit_mat[,j] - rna_common2[,j])^2)
  flat_mse <- sd(rna_common2[,j])
  abs(unimodal_mse - flat_mse)/unimodal_mse
})
atac_quality <- sapply(1:p, function(j){
  unimodal_mse <- mean((atac_unimodal_fit_mat[,j] - atac_everything2[,j])^2)
  flat_mse <- sd(atac_everything2[,j])
  abs(unimodal_mse - flat_mse)/unimodal_mse
})

num_considered <- 100
idx2 <- order(rna_quality, decreasing = T)[1:num_considered]
idx3 <- order(atac_quality, decreasing = T)[1:num_considered]
idx <- intersect(idx2, idx3)
length(idx)
col_vec <- grDevices::colorRampPalette(c('lightgrey', 'blue'))(nrow(rna_common2))
png("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_activity_everythingATAC2commonRNA_individualgenes_unimodal_downward.png",
    height = 2000, width = 3000, units = "px", res = 300)
par(mfrow = c(3,5), mar = c(4,4,4,0.5))
for(j in idx[1:15]){
  shuff_idx <- sample(1:n)
  plot(atac_everything2[shuff_idx,j], rna_common2[shuff_idx,j],
       main = paste0("Correlation: ", round(cor_vec[j], 2)),
       xlab = paste0("Gene activity for ", colnames(rna_common)[j]),
       ylab = paste0("Common RNA for ", colnames(rna_common2)[j]),
       col = col_vec[shuff_idx], pch = 16, cex = 0.5
  )
}
graphics.off()
png("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_activity_everythingATAC2commonRNA_individualgenes_unimodal_downward2.png",
    height = 2000, width = 3000, units = "px", res = 300)
par(mfrow = c(3,5), mar = c(4,4,4,0.5))
for(j in idx[1:15]){
  plot(NA, 
       xlim = c(1,n),
       ylim = range(c(atac_everything2[,j], rna_common2[,j])),
       xlab = "Rank (via Slingshot)",
       ylab = paste0("Expression for ", colnames(rna_common2)[j])
  )
  
  x_vec <- c(1:n, 1:n)
  y_vec <- c(rna_common2[,j], atac_everything2[,j])
  col_vec <- c(rep(rgb(0.5, 0.5, 0.5, 0.5), n),
               rep(rgb(0.75, 0, 0, 0.5), n))
  shuff_idx <- sample(1:length(x_vec))
  points(x = x_vec[shuff_idx], y = y_vec[shuff_idx], 
         col = col_vec[shuff_idx], pch = 16, cex = 0.5)
  lines(x = 1:n, y = rna_unimodal_fit_mat[,j],
        col = "white", lwd = 5)
  lines(x = 1:n, y = rna_unimodal_fit_mat[,j],
        col = "black", lwd = 3)
  lines(x = 1:n, y = atac_unimodal_fit_mat[,j],
        col = "white", lwd = 5)
  lines(x = 1:n, y = atac_unimodal_fit_mat[,j],
        col = 2, lwd = 3)
}
graphics.off()


col_pal <- grDevices::colorRampPalette(c('lightgrey', 'blue'))(100)
seq_vec <- seq(min(alignment_vec), max(alignment_vec), length.out = 100)
col_vec <- sapply(alignment_vec, function(val){
  col_pal[which.min(abs(val - seq_vec))]
})
png("../../../../out/figures/Writeup14m/Writeup14m_10x_greenleaf_activity_everythingATAC2commonRNA_individualgenes_unimodal_downward_alignment.png",
    height = 2000, width = 3000, units = "px", res = 300)
par(mfrow = c(3,5), mar = c(4,4,4,0.5))
for(j in idx[1:15]){
  plot(atac_everything2[,j], rna_common2[,j],
       xlab = paste0("Gene activity for ", colnames(rna_common)[j]),
       ylab = paste0("Common RNA for ", colnames(rna_common2)[j]),
       col = col_vec, pch = 16, cex = 0.5
  )
}
graphics.off()
