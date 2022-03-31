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

####

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

####

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

####

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

####

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

####

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
gene_idx_upward <- intersect(idx2, idx3)
length(gene_idx_upward)

####

p <- ncol(rna_common2); n <- nrow(rna_common2)
rna_decreasing_fit_mat <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  UniIsoRegression::reg_1d(rna_common2[,j],
                           w_vec = rep(1, n),
                           metric = 2,
                           unimodal = F,
                           decreasing = T)
})
atac_decreasing_fit_mat <- sapply(1:p, function(j){
  if(j %% floor(p/10) == 0) cat('*')
  
  UniIsoRegression::reg_1d(atac_everything2[,j],
                           w_vec = rep(1, n),
                           metric = 2,
                           unimodal = F,
                           decreasing = T)
})
rna_quality <- sapply(1:p, function(j){
  decreasing_mse <- mean((rna_decreasing_fit_mat[,j] - rna_common2[,j])^2)
  flat_mse <- sd(rna_common2[,j])
  abs(decreasing_mse - flat_mse)/decreasing_mse
})
atac_quality <- sapply(1:p, function(j){
  decreasing_mse <- mean((atac_decreasing_fit_mat[,j] - atac_everything2[,j])^2)
  flat_mse <- sd(atac_everything2[,j])
  abs(decreasing_mse - flat_mse)/decreasing_mse
})

num_considered <- 100
idx2 <- order(rna_quality, decreasing = T)[1:num_considered]
idx3 <- order(atac_quality, decreasing = T)[1:num_considered]
gene_idx_downward <- intersect(idx2, idx3)
length(gene_idx_downward)

save(rna_decreasing_fit_mat, atac_decreasing_fit_mat,
     rna_unimodal_fit_mat, atac_unimodal_fit_mat,
     gene_idx_downward, gene_idx_upward,
     rna_common2, atac_everything2,
     file = "../../../../out/Writeup14m/Writeup14m_10x_greenleaf_geneactivity_grangerTest.RData")

####

## now try lmtest::grangertest
granger_vec <- sapply(gene_idx, function(j){
  df <- data.frame(vec1 = atac_unimodal_fit_mat[,j],
                   vec2 = rna_unimodal_fit_mat[,j])
 
  tmp <- lmtest::grangertest(vec2 ~ vec1,
                             order = 1, 
                             data = df)
  tmp["Pr(>F)"][2,1]
})
round(granger_vec, 2)

granger_vec <- sapply(gene_idx, function(j){
  df <- data.frame(vec1 = atac_unimodal_fit_mat[,j],
                   vec2 = rna_unimodal_fit_mat[,j])
  
  tmp <- lmtest::grangertest(vec1 ~ vec2,
                             order = 1, 
                             data = df)
  tmp["Pr(>F)"][2,1]
})
round(granger_vec, 2)

###########################

