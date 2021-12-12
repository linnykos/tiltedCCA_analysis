rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_bm25_preprocessed.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(bm)

mat_1 <- t(bm[["RNA"]]@scale.data)
mat_2 <- t(bm[["ADT"]]@scale.data)

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

#########################

rank_1 <- 30; rank_2 <- 18; nn <- 30
set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1b, mat_2b, 
                                      dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      num_neigh = nn, 
                                      metacell_clustering_1 = NA,
                                      metacell_clustering_2 = NA,
                                      fix_tilt_perc = 0.5, verbose = T)

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)

n <- nrow(bm@meta.data)
svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_1, 
                                       K = rank_1, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_2, K = rank_2, K_full_rank = F,
                                       rescale = F,
                                       mean_vec = F, sd_vec = F,
                                       symmetric = F)
dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
l2_vec <- apply(dimred_1, 1, function(x){multiomicCCA:::.l2norm(x)})
dimred_1 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_1)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
l2_vec <- apply(dimred_2, 1, function(x){multiomicCCA:::.l2norm(x)})
dimred_2 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_2)

set.seed(10)
common_umap <- Seurat::RunUMAP(cbind(dimred_1, dimred_2), 
                               metric = "euclidean",
                               reduction.key = "umapCommon_")
rownames(common_umap@cell.embeddings) <- rownames(bm@meta.data)
bm[["dcca_common"]] <- Seurat::CreateDimReducObject(common_umap@cell.embeddings)

######################

dims_1 <- 1:30; dims_2 <- 1:18
rank_1 <- max(dims_1); rank_2 <- max(dims_2)
n <- nrow(mat_1)

svd_1 <- multiomicCCA:::.svd_truncated(mat_1, K = rank_1, symmetric = F, rescale = F, 
                                       mean_vec = F, sd_vec = F, K_full_rank = F)
svd_2 <- multiomicCCA:::.svd_truncated(mat_2, K = rank_2, symmetric = F, rescale = F, 
                                       mean_vec = F, sd_vec = F, K_full_rank = F)

svd_1 <- multiomicCCA:::.check_svd(svd_1, dims = dims_1)
svd_2 <- multiomicCCA:::.check_svd(svd_2, dims = dims_2)

dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)

r <- 18
basis_list <- lapply(1:r, function(k){
  multiomicCCA:::.representation_2d(dcca_res$score_1[,k], dcca_res$score_2[,k])
})

circle_list <- lapply(1:r, function(k){
  vec1 <- basis_list[[k]]$rep1
  vec2 <- basis_list[[k]]$rep2
  multiomicCCA:::.construct_circle(vec1, vec2)
})

radian_vec <- sapply(1:r, function(k){
  multiomicCCA:::.compute_radian(circle = circle_list[[k]],
                                 enforce_boundary = T,
                                 percentage_val = 0.5, 
                                 vec1 = basis_list[[k]]$rep1,
                                 vec2 = basis_list[[k]]$rep2)
})

common_representation <- sapply(1:r, function(k){
  multiomicCCA:::.position_from_circle(circle_list[[k]], radian_vec[k])
})

common_score <- sapply(1:r, function(k){
  basis_list[[k]]$basis_mat %*% common_representation[,k]
})

common_mat <-  multiomicCCA:::.convert_common_score_to_mat(common_score,
                                                           dcca_res$score_1,
                                                           dcca_res$score_2,
                                                           dcca_res$svd_1, 
                                                           dcca_res$svd_2)


candidate_vec <- which(bm@meta.data[,"celltype.l2"] == "CD4 Naive")
rann_res_1 <- RANN::nn2(data = dimred_1, query = dimred_1[candidate_vec,], 
                        k = 31)$nn.idx[,-1]
rann_res_2 <- RANN::nn2(data = dimred_2, query = dimred_2[candidate_vec,], 
                        k = 31)$nn.idx[,-1]
rann_res_3 <- RANN::nn2(data = common_mat, query = common_mat[candidate_vec,], 
                        k = 31)$nn.idx[,-1]
# remove entries where ADT has CD8 naive neighbors
bool_vec <- sapply(1:length(candidate_vec), function(i){
  any(bm$celltype.l2[rann_res_2[i,]] == "CD8 Naive")
})
if(any(bool_vec)){
  candidate_vec <- candidate_vec[!bool_vec]
  rann_res_1 <- rann_res_1[!bool_vec,]
  rann_res_2 <- rann_res_2[!bool_vec,]
  rann_res_3 <- rann_res_3[!bool_vec,]
}

overlap_vec <- sapply(1:length(candidate_vec), function(i){
  min(c(length(intersect(rann_res_1[i,], rann_res_3[i,])),
        length(intersect(rann_res_2[i,], rann_res_3[i,]))))
})
max(overlap_vec)
target_idx <- which.max(overlap_vec)
cell_idx <- candidate_vec[target_idx]
table(bm@meta.data[rann_res_1[target_idx,],"celltype.l2"])
table(bm@meta.data[rann_res_2[target_idx,],"celltype.l2"])
length(intersect(rann_res_1[target_idx,], rann_res_3[target_idx,]))
length(intersect(rann_res_2[target_idx,], rann_res_3[target_idx,]))

##############################3

anchor_name <- "rna.umap"
other_names <- c("adt.umap", "dcca_common")

for(umap_name in other_names){
  print(umap_name)
  u_mat1 <- bm[[anchor_name]]@cell.embeddings
  u_mat2 <- bm[[umap_name]]@cell.embeddings
  tmp <- svd(t(u_mat1) %*% u_mat2)
  rotation_mat <- tmp$u %*% t(tmp$v)
  tmp <- u_mat2 %*% t(rotation_mat)
  rownames(tmp) <- rownames(bm@meta.data)
  colnames(tmp) <- colnames(bm[[umap_name]]@cell.embeddings)
  bm[[umap_name]]@cell.embeddings <- tmp
}

############

png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_rna_neigh.png",
    height = 2050, width = 1800, res = 300, units = "px")
tmp <- bm[["rna.umap"]]@cell.embeddings
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_1[target_idx,],1], tmp[rann_res_1[target_idx,],2], col = "white", pch = 16, cex = 4)
points(tmp[rann_res_1[target_idx,],1], tmp[rann_res_1[target_idx,],2], col = 2, pch = 16, cex = 3)
graphics.off()

png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_adt_neigh.png",
    height = 2050, width = 1800, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
tmp <- bm[["adt.umap"]]@cell.embeddings
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_2[target_idx,],1], tmp[rann_res_2[target_idx,],2], col = "white", pch = 16, cex = 4)
points(tmp[rann_res_2[target_idx,],1], tmp[rann_res_2[target_idx,],2], col = 2, pch = 16, cex = 3)
graphics.off()

##########

color_vec3 <- c(rgb(127, 127, 127, max = 255), 
                rgb(212, 65, 86, max = 255), 
                rgb(72, 161, 217, max = 255), 
                rgb(249, 97, 255, max = 255), 
                rgb(77, 194, 61, max = 255))
names(color_vec3) <- c("none", "neigh", "neighdcca", "intersect", "target")


png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_rna_onlycell.png",
    height = 2050, width = 1800, res = 300, units = "px")
tmp <- bm[["rna.umap"]]@cell.embeddings
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()

png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_adt_onlycell.png",
    height = 2050, width = 1800, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
tmp <- bm[["adt.umap"]]@cell.embeddings
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()

png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_dcca_neighandcell.png",
    height = 2050, width = 1800, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
tmp <- bm[["dcca_common"]]@cell.embeddings
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = "white", pch = 16, cex = 4)
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = color_vec3["neighdcca"], pch = 16, cex = 3)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()

png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_rna_neighandcell.png",
    height = 2050, width = 1800, res = 300, units = "px")
tmp <- bm[["rna.umap"]]@cell.embeddings
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_1[target_idx,],1], tmp[rann_res_1[target_idx,],2], col = "white", pch = 16, cex = 4)
points(tmp[rann_res_1[target_idx,],1], tmp[rann_res_1[target_idx,],2], col = color_vec3["neigh"], pch = 16, cex = 3)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()

png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_adt_neighandcell.png",
    height = 2050, width = 1800, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
tmp <- bm[["adt.umap"]]@cell.embeddings
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_2[target_idx,],1], tmp[rann_res_2[target_idx,],2], col = "white", pch = 16, cex = 4)
points(tmp[rann_res_2[target_idx,],1], tmp[rann_res_2[target_idx,],2], col = color_vec3["neigh"], pch = 16, cex = 3)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()

##############

png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_rna_neighintersect.png",
    height = 2050, width = 1800, res = 300, units = "px")
tmp <- bm[["rna.umap"]]@cell.embeddings
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_1[target_idx,],1], tmp[rann_res_1[target_idx,],2], col = "gray", pch = 16, cex = 4)
points(tmp[rann_res_1[target_idx,],1], tmp[rann_res_1[target_idx,],2], col = color_vec3["neigh"], pch = 16, cex = 3)
tmp_idx <- intersect(rann_res_1[target_idx,], rann_res_3[target_idx,])
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = color_vec3["intersect"], pch = 16, cex = 3)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()

png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_adt_neighintersect.png",
    height = 2050, width = 1800, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
tmp <- bm[["adt.umap"]]@cell.embeddings
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_2[target_idx,],1], tmp[rann_res_2[target_idx,],2], col = "gray", pch = 16, cex = 4)
points(tmp[rann_res_2[target_idx,],1], tmp[rann_res_2[target_idx,],2], col = color_vec3["neigh"], pch = 16, cex = 3)
tmp_idx <- intersect(rann_res_2[target_idx,], rann_res_3[target_idx,])
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = color_vec3["intersect"], pch = 16, cex = 3)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()


png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_dcca_neighintersect1.png",
    height = 2050, width = 1800, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
tmp <- bm[["dcca_common"]]@cell.embeddings
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = "gray", pch = 16, cex = 4)
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = color_vec3["neighdcca"], pch = 16, cex = 3)
tmp_idx <- intersect(rann_res_1[target_idx,], rann_res_3[target_idx,])
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = color_vec3["intersect"], pch = 16, cex = 3)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()


png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_dcca_neighintersect2.png",
    height = 2050, width = 1800, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
tmp <- bm[["dcca_common"]]@cell.embeddings
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = "gray", pch = 16, cex = 4)
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = color_vec3["neighdcca"], pch = 16, cex = 3)
tmp_idx <- intersect(rann_res_2[target_idx,], rann_res_3[target_idx,])
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = color_vec3["intersect"], pch = 16, cex = 3)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()

#############################################

tmp <- bm[["adt.umap"]]@cell.embeddings
zz <- tmp[c(rann_res_2[target_idx,],cell_idx),]
xlim <- range(zz[,1]); ylim <- range(zz[,2])
png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_adt_neighintersect_zoomin.png",
    height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xlim = xlim, ylim = ylim, asp = T,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_2[target_idx,],1], tmp[rann_res_2[target_idx,],2], col = "gray", pch = 16, cex = 5)
points(tmp[rann_res_2[target_idx,],1], tmp[rann_res_2[target_idx,],2], col = color_vec3["neigh"], pch = 16, cex = 4)
tmp_idx <- intersect(rann_res_2[target_idx,], rann_res_3[target_idx,])
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = "white", pch = 16, cex = 5)
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = color_vec3["intersect"], pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 5)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 4)
graphics.off()


tmp <- bm[["dcca_common"]]@cell.embeddings
zz <- tmp[c(rann_res_3[target_idx,],cell_idx),]
xlim <- range(zz[,1]); ylim <- range(zz[,2])
png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_dcca_neighintersect1_zoomin1.png",
    height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xlim = xlim, ylim = ylim, asp = T,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = "white", pch = 16, cex = 4)
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = color_vec3["neighdcca"], pch = 16, cex = 3)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()

tmp <- bm[["dcca_common"]]@cell.embeddings
zz <- tmp[c(rann_res_3[target_idx,],cell_idx),]
xlim <- range(zz[,1]); ylim <- range(zz[,2])
png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_dcca_neighintersect1_zoomin2.png",
    height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xlim = xlim, ylim = ylim, asp = T,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = "gray", pch = 16, cex = 4)
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = color_vec3["neighdcca"], pch = 16, cex = 3)
tmp_idx <- intersect(rann_res_1[target_idx,], rann_res_3[target_idx,])
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = color_vec3["intersect"], pch = 16, cex = 3)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()


tmp <- bm[["dcca_common"]]@cell.embeddings
zz <- tmp[c(rann_res_3[target_idx,],cell_idx),]
xlim <- range(zz[,1]); ylim <- range(zz[,2])
png("../../../../out/figures/Writeup14g/Writeup14g_citeseq_bm_dcca_neighintersect2_zoomin2.png",
    height = 1000, width = 1000, res = 300, units = "px")
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(tmp[,1], tmp[,2], 
     pch = 16, col = "gray50", cex = 0.5,
     xlim = xlim, ylim = ylim, asp = T,
     xaxt = "n", yaxt = "n", bty = "n")
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = "gray", pch = 16, cex = 4)
points(tmp[rann_res_3[target_idx,],1], tmp[rann_res_3[target_idx,],2], col = color_vec3["neighdcca"], pch = 16, cex = 3)
tmp_idx <- intersect(rann_res_2[target_idx,], rann_res_3[target_idx,])
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[tmp_idx,1], tmp[tmp_idx,2], col = color_vec3["intersect"], pch = 16, cex = 3)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = "white", pch = 16, cex = 4)
points(tmp[cell_idx,1], tmp[cell_idx,2], col = color_vec3["target"], pch = 16, cex = 3)
graphics.off()

