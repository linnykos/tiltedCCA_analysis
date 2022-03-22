rm(list=ls())
load("../../../../out/Writeup14l/Writeup14l_10x_greenleaf2_tcca.RData")
library(Seurat)

set.seed(10)
greenleaf <- Seurat::FindMultiModalNeighbors(greenleaf, reduction.list = list("pca", "lsi"), 
                                             dims.list = list(1:50, 2:50))
set.seed(10)
greenleaf <- Seurat::FindClusters(greenleaf,
                                  graph.name = "wsnn", algorithm = 3, 
                                  resolution = 2)

plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         group.by = "seurat_clusters", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nSeurat WNN clusters, on Common embedding"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_seuratWNNcluster.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

cell_names <- colnames(greenleaf)[which(greenleaf$seurat_clusters %in% as.character(c(13,15,1,0,3,31,19,2,11,18,12,8,28)))]
plot1 <- Seurat::DimPlot(greenleaf, reduction = "common_tcca",
                         cells.highlight = cell_names)
plot1 <- plot1 + ggplot2::ggtitle(paste0("Human brain (10x, RNA+ATAC)\nSeurat WNN clusters (subset)"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_seuratWNNcluster_subset.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

#########################

multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
dimred_1 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, 
                                     apply_postDimred = T, what = "common_mat")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred_2 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, what = "common_dimred")
dimred <- cbind(dimred_1, dimred_2)
svd_tmp <- irlba::irlba(dimred, nv = 25)
dimred <- tiltedCCA:::.mult_mat_vec(svd_tmp$u[,1:10], svd_tmp$d[1:10])
lineage_order <- as.character(c(13,15,1,0,3,31,19,2,11,18,12,8,28))
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
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_pseudotime.png"),
                plot3, device = "png", width = 10, height = 5, units = "in")


png("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_pseudotime_values.png",
    height = 1500, width = 1500, units = "px", res = 300)
plot(sort(pseudotime_vec), xlab = "Index (sorted order)", ylab = "Pseudotime",
     main = "Pseudotime (ordered) of GluN4 lineage",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
graphics.off()

###############################

multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
dimred_1 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, 
                                     what = "common_mat")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
dimred_2 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = T, 
                                     what = "common_dimred")
dimred <- cbind(dimred_1, dimred_2)
snn_mat <- tiltedCCA:::.form_snn_mat(mat = dimred,
                                     num_neigh = multiSVD_obj$param$snn_num_neigh,
                                     bool_cosine = multiSVD_obj$param$snn_bool_cosine, 
                                     bool_intersect = multiSVD_obj$param$snn_bool_intersect, 
                                     min_deg = multiSVD_obj$param$snn_min_deg,
                                     tol = 1e-4,
                                     verbose = 0)

###############################

rna_common <- multiSVD_obj$common_mat_1
n <- nrow(rna_common)
Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
mat_1 <- scale(mat_1)

fate_res <- .compute_fate_direction(mat_firsttime = rna_common,
                                    mat_secondtime = mat_1,
                                    pseudotime_rank = greenleaf$pseudotime_rank,
                                    snn_mat = snn_mat,
                                    verbose = 1)
fate_vec_full <- fate_res$fate_vec
idx <- which(!is.na(fate_vec_full))
fate_vec_full[idx] <- sign(fate_vec_full[idx])*abs(fate_vec_full[idx])^(1/4)
col_spacing <- c(seq(-1,0,length.out=21),seq(0,1,length.out=21)[-1])
col_palette <- c(grDevices::colorRampPalette(c(2, "white"))(21),
                 grDevices::colorRampPalette(c("white", 3))(21)[-1])
col_vec <- sapply(fate_vec_full, function(val){
  if(is.na(val)) return(NA)
  col_palette[which.min(abs(val - col_spacing))]
})

png("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_fateDirection.png",
    height = 3000, width = 3000, units = "px", res = 300)
umap_mat <- greenleaf[["common_tcca"]]@cell.embeddings
plot(umap_mat[,1], umap_mat[,2], asp = T,
     xlab = "common_tcca_1", ylab = "common_tcca_2",
     main = "Fate direction (Nearby RNA expression\nregressed onto cell's common RNA)",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
idx <- which(!is.na(col_vec))
points(umap_mat[idx,1], umap_mat[idx,2],
       pch = 16, col = col_vec[idx], cex = 0.5)
graphics.off()



#########################################

pseudotime_rank <- greenleaf$pseudotime_rank
snn_1 <- multiSVD_obj$snn_list$snn_1
snn_2 <- multiSVD_obj$snn_list$snn_2
cell_idx <- which(!is.na(pseudotime_rank))

fate_vec <- sapply(1:length(cell_idx), function(k){
  if(k %% floor(length(cell_idx)/10) == 0) cat('*')
  
  i <- cell_idx[k]
  nn_idx_2 <- tiltedCCA:::.nonzero_col(mat = snn_2, col_idx = i, bool_value = F)
  nn_idx_2 <- intersect(nn_idx_2, cell_idx)
  if(length(nn_idx_2) == 0) return(NA)
  cell_rank <- pseudotime_rank[i]
  nn_rank <- pseudotime_rank[nn_idx_2]
  nn_val <- nn_rank-cell_rank
  nn_val <- nn_val/max(abs(nn_val))
  
  jaccard_vec <- sapply(nn_idx_2, function(j){
    nn_idx_1 <- tiltedCCA:::.nonzero_col(mat = snn_1, col_idx = j, bool_value = F)
    nn_idx_1 <- intersect(nn_idx_1, cell_idx)
    
    length(intersect(nn_idx_2, nn_idx_1))/length(unique(c(nn_idx_2, nn_idx_1)))
  })
  nn_val[which.max(jaccard_vec)]
})

fate_vec_full <- rep(NA, n)
fate_vec_full[which(!is.na(pseudotime_rank))] <- fate_vec
quantile(fate_vec_full, probs = seq(0,1,length.out=11), na.rm = T)
quantile(fate_vec_full, probs = c(0.01,.99), na.rm = T)
idx <- which(!is.na(fate_vec_full))
fate_vec_full[idx] <- sign(fate_vec_full[idx])*abs(fate_vec_full[idx])^(1/4)
col_spacing <- c(seq(-1,0,length.out=21),seq(0,1,length.out=21)[-1])
col_palette <- c(grDevices::colorRampPalette(c(2, "white"))(21),
                 grDevices::colorRampPalette(c("white", 3))(21)[-1])
col_vec <- sapply(fate_vec_full, function(val){
  if(is.na(val)) return(NA)
  col_palette[which.min(abs(val - col_spacing))]
})

png("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_fateDirection_jaccard.png",
    height = 3000, width = 3000, units = "px", res = 300)
umap_mat <- greenleaf[["common_tcca"]]@cell.embeddings
plot(umap_mat[,1], umap_mat[,2], asp = T,
     xlab = "common_tcca_1", ylab = "common_tcca_2",
     main = "Fate direction (Nearby RNA NN's\nwith highest overlap with ATAC's NN)",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
idx <- which(!is.na(col_vec))
points(umap_mat[idx,1], umap_mat[idx,2],
       pch = 16, col = col_vec[idx], cex = 0.5)
graphics.off()

####################################

# let's try something else
Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
mat_1 <- scale(mat_1)
common_score <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, 
                                         what = "common_score")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
score_1 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F,
                                    what = "score")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
score_2 <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, 
                                    what = "score")
k <- ncol(common_score)
n <- nrow(common_score)
# rna_dimred <- common_score %*% crossprod(score_1[,1:k], tiltedCCA:::.mult_mat_vec(svd_1$u, svd_1$d))/n
# rna_dimred <- tiltedCCA:::.mult_mat_vec(svd_1$u, svd_1$d)/n
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
# atac_dimred <- common_score %*% crossprod(score_2[,1:k], tiltedCCA:::.mult_mat_vec(svd_2$u, svd_2$d)) %*% rotation_mat/n
tmp <- common_score %*% crossprod(score_2[,1:k], svd_2$u) %*% rotation_mat/n
atac_pred <-  tcrossprod(tiltedCCA:::.mult_mat_vec(tmp, svd_1$d), svd_1$v)

fate_res <- .compute_fate_direction(mat_firsttime = atac_pred,
                                    mat_secondtime = mat_1,
                                    pseudotime_rank = greenleaf$pseudotime_rank,
                                    snn_mat = snn_mat,
                                    verbose = 1)
fate_vec_full <- fate_res$fate_vec
idx <- which(!is.na(fate_vec_full))
fate_vec_full[idx] <- sign(fate_vec_full[idx])*abs(fate_vec_full[idx])^(1/4)
col_spacing <- c(seq(-1,0,length.out=21),seq(0,1,length.out=21)[-1])
col_palette <- c(grDevices::colorRampPalette(c(2, "white"))(21),
                 grDevices::colorRampPalette(c("white", 3))(21)[-1])
col_vec <- sapply(fate_vec_full, function(val){
  if(is.na(val)) return(NA)
  col_palette[which.min(abs(val - col_spacing))]
})

png("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_fateDirection_commonAtac2everythingRna.png",
    height = 3000, width = 3000, units = "px", res = 300)
umap_mat <- greenleaf[["common_tcca"]]@cell.embeddings
plot(umap_mat[,1], umap_mat[,2], asp = T,
     xlab = "common_tcca_1", ylab = "common_tcca_2",
     main = "Fate direction (Nearby everything RNA's\npredicted by common ATAC)",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
idx <- which(!is.na(col_vec))
points(umap_mat[idx,1], umap_mat[idx,2],
       pch = 16, col = col_vec[idx], cex = 0.5)
graphics.off()

#########################################

# let's do the analysis with only the important genes
rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
mat_1 <- rna_common + rna_distinct
p <- ncol(rna_common)
gene_alignment <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  df <- cbind(mat_1[,j], multiSVD_obj$laplacian_list$common_laplacian)
  df <- as.data.frame(df)
  colnames(df)[1] <- "y"
  lm_res <- stats::lm(y ~ ., data = df)
  summary(lm_res)$r.squared
})
gene_var <- matrixStats::colSds(mat_1)
names(gene_var) <- colnames(mat_1)

gene_names <- intersect(colnames(mat_1)[which(gene_alignment >= 0.6)],
                        colnames(mat_1)[which(gene_var >= 0.4)])

fate_vec <- sapply(1:length(cell_idx), function(k){
  if(k %% floor(length(cell_idx)/10) == 0) cat('*')
  
  i <- cell_idx[k]
  nn_idx <- tiltedCCA:::.nonzero_col(mat = snn_1, col_idx = i, bool_value = F)
  nn_idx <- intersect(nn_idx, cell_idx)
  if(length(nn_idx) == 0) return(NA)
  cell_rank <- pseudotime_rank[i]
  nn_rank <- pseudotime_rank[nn_idx]
  nn_val <- nn_rank-cell_rank
  nn_val <- nn_val/max(abs(nn_val))
  
  fit_vec <- sapply(nn_idx, function(j){
    df <- data.frame(rna = mat_1[j,gene_names],
                     atac = atac_pred[i,gene_names])
    lm_res <- stats::lm(rna ~ atac, data = df)
    summary(lm_res)$r.squared
  })
  nn_val[which.max(fit_vec)]
})

fate_vec_full <- rep(NA, n)
fate_vec_full[which(!is.na(pseudotime_rank))] <- fate_vec
quantile(fate_vec_full, probs = seq(0,1,length.out=11), na.rm = T)
quantile(fate_vec_full, probs = c(0.01,.99), na.rm = T)
fate_vec_full[idx] <- sign(fate_vec_full[idx])*abs(fate_vec_full[idx])^(1/4)
idx <- which(!is.na(fate_vec_full))
col_spacing <- c(seq(-1,0,length.out=21),seq(0,1,length.out=21)[-1])
col_palette <- c(grDevices::colorRampPalette(c(2, "white"))(21),
                 grDevices::colorRampPalette(c("white", 3))(21)[-1])
col_vec <- sapply(fate_vec_full, function(val){
  if(is.na(val)) return(NA)
  col_palette[which.min(abs(val - col_spacing))]
})

png("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_fateDirection_commonAtac2everythingRna_genesubset.png",
    height = 3000, width = 3000, units = "px", res = 300)
umap_mat <- greenleaf[["common_tcca"]]@cell.embeddings
plot(umap_mat[,1], umap_mat[,2], asp = T,
     xlab = "common_tcca_1", ylab = "common_tcca_2",
     main = "Fate direction (Nearby everything RNA's\npredicted by common ATAC)",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
idx <- which(!is.na(col_vec))
points(umap_mat[idx,1], umap_mat[idx,2],
       pch = 16, col = col_vec[idx], cex = 0.5)
graphics.off()

###################################

Seurat::DefaultAssay(greenleaf) <- "SCT"
mat_1 <- Matrix::t(greenleaf[["SCT"]]@data[Seurat::VariableFeatures(object = greenleaf),])
mat_1 <- scale(mat_1)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)

fate_res <- .compute_fate_direction(mat_firsttime = atac_pred,
                                    mat_secondtime = mat_1,
                                    pseudotime_rank = greenleaf$pseudotime_rank,
                                    snn_mat = snn_mat,
                                    verbose = 1)
fate_vec_full <- fate_res$fate_vec
idx <- which(!is.na(fate_vec_full))
fate_vec_full[idx] <- sign(fate_vec_full[idx])*abs(fate_vec_full[idx])^(1/4)
col_spacing <- c(seq(-1,0,length.out=21),seq(0,1,length.out=21)[-1])
col_palette <- c(grDevices::colorRampPalette(c(2, "white"))(21),
                 grDevices::colorRampPalette(c("white", 3))(21)[-1])
col_vec <- sapply(fate_vec_full, function(val){
  if(is.na(val)) return(NA)
  col_palette[which.min(abs(val - col_spacing))]
})

png("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_fateDirection_everythingAtac2everythingRna.png",
    height = 3000, width = 3000, units = "px", res = 300)
umap_mat <- greenleaf[["common_tcca"]]@cell.embeddings
plot(umap_mat[,1], umap_mat[,2], asp = T,
     xlab = "common_tcca_1", ylab = "common_tcca_2",
     main = "Fate direction (Nearby everything RNA's\npredicted by everything ATAC)",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
idx <- which(!is.na(col_vec))
points(umap_mat[idx,1], umap_mat[idx,2],
       pch = 16, col = col_vec[idx], cex = 0.5)
graphics.off()

#################################

rna_common <- multiSVD_obj$common_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_pred <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)

fate_res <- .compute_fate_direction(mat_firsttime = atac_pred,
                                    mat_secondtime = rna_common,
                                    pseudotime_rank = greenleaf$pseudotime_rank,
                                    snn_mat = snn_mat,
                                    verbose = 1)
fate_vec_full <- fate_res$fate_vec
idx <- which(!is.na(fate_vec_full))
fate_vec_full[idx] <- sign(fate_vec_full[idx])*abs(fate_vec_full[idx])^(1/4)
idx <- which(!is.na(fate_vec_full))
col_spacing <- c(seq(-1,0,length.out=21),seq(0,1,length.out=21)[-1])
col_palette <- c(grDevices::colorRampPalette(c(2, "white"))(21),
                 grDevices::colorRampPalette(c("white", 3))(21)[-1])
col_vec <- sapply(fate_vec_full, function(val){
  if(is.na(val)) return(NA)
  col_palette[which.min(abs(val - col_spacing))]
})

png("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_fateDirection_everythingAtac2commonRna.png",
    height = 3000, width = 3000, units = "px", res = 300)
umap_mat <- greenleaf[["common_tcca"]]@cell.embeddings
plot(umap_mat[,1], umap_mat[,2], asp = T,
     xlab = "common_tcca_1", ylab = "common_tcca_2",
     main = "Fate direction (Nearby common RNA's\npredicted by everything ATAC)",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
idx <- which(!is.na(col_vec))
points(umap_mat[idx,1], umap_mat[idx,2],
       pch = 16, col = col_vec[idx], cex = 0.5)
graphics.off()


fate_vec_full <- rep(NA, n)
fate_vec_full[which(!is.na(pseudotime_rank))] <- fate_vec
quantile(fate_vec_full, probs = seq(0,1,length.out=11), na.rm = T)
quantile(fate_vec_full, probs = c(0.01,.99), na.rm = T)
fate_vec_full[idx] <- sign(fate_vec_full[idx])*abs(fate_vec_full[idx])^(1/2)
idx <- which(!is.na(fate_vec_full))
col_spacing <- c(seq(-1,0,length.out=21),seq(0,1,length.out=21)[-1])
col_palette <- c(grDevices::colorRampPalette(c(2, "white"))(21),
                 grDevices::colorRampPalette(c("white", 3))(21)[-1])
col_vec <- sapply(fate_vec_full, function(val){
  if(is.na(val)) return(NA)
  col_palette[which.min(abs(val - col_spacing))]
})

png("../../../../out/figures/Writeup14l/Writeup14l_10x_greenleaf2_fateDirection_everythingAtac2commonRna2.png",
    height = 3000, width = 3000, units = "px", res = 300)
umap_mat <- greenleaf[["common_tcca"]]@cell.embeddings
plot(umap_mat[,1], umap_mat[,2], asp = T,
     xlab = "common_tcca_1", ylab = "common_tcca_2",
     main = "Fate direction (Nearby common RNA's\npredicted by everything ATAC)",
     pch = 16, col = rgb(0.5,0.5,0.5,0.5))
idx <- which(!is.na(col_vec))
points(umap_mat[idx,1], umap_mat[idx,2],
       pch = 16, col = col_vec[idx], cex = 0.5)
graphics.off()

