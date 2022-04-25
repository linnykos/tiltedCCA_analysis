rm(list=ls())
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")

library(Seurat); library(Signac)

input_obj <- multiSVD_obj

input_obj <- tiltedCCA:::.set_defaultAssay(input_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(input_obj)
input_obj <- tiltedCCA:::.set_defaultAssay(input_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(input_obj)

n <- nrow(svd_1$u)
metacell_clustering_list <- tiltedCCA:::.get_metacell(input_obj,
                                                      resolution = "cell", 
                                                      type = "list", 
                                                      what = "metacell_clustering")
averaging_mat <- tiltedCCA:::.generate_averaging_matrix(metacell_clustering_list = metacell_clustering_list,
                                                        n = n)

target_dimred <- tiltedCCA:::.get_Laplacian(input_obj, bool_common = T)
param <- tiltedCCA:::.get_param(input_obj)
snn_bool_cosine <- param$snn_bool_cosine
snn_bool_intersect <- param$snn_bool_intersect
snn_k <- param$snn_latent_k
snn_min_deg <- param$snn_min_deg
snn_num_neigh <- param$snn_num_neigh

cca_res <- tiltedCCA:::.cca(svd_1, svd_2, dims_1 = NA, dims_2 = NA, return_scores = F)

###########

full_rank <- length(cca_res$obj_vec)
n <- nrow(svd_1$u)
tmp <- tiltedCCA:::.compute_unnormalized_scores(svd_1, svd_2, cca_res)
score_1 <- tmp$score_1; score_2 <- tmp$score_2
stopifnot(ncol(score_1) == length(svd_1$d), ncol(score_2) == length(svd_2$d),
          nrow(score_1) == nrow(score_2))
obj_vec <- cca_res$obj_vec

##########

rank_c <- min(ncol(score_1), ncol(score_2))
stopifnot(all(sapply(1:rank_c, function(k){
  val <- score_1[,k] %*% score_2[,k]; val >= 0 
}))) # ensures score matrices contain pair of acute vectors

basis_list <- lapply(1:rank_c, function(k){
  tiltedCCA:::.representation_2d(score_1[,k], score_2[,k])
})

circle_list <- lapply(1:rank_c, function(k){
  vec1 <- basis_list[[k]]$rep1
  vec2 <- basis_list[[k]]$rep2
  tiltedCCA:::.construct_circle(vec1, vec2)
})

###############

r <- length(basis_list)
enforce_boundary <- F

percentage <- multiSVD_obj$tcca_obj$tilt_perc
radian_vec <- sapply(1:r, function(k){
  tiltedCCA:::.compute_radian(circle = circle_list[[k]],
                              enforce_boundary = enforce_boundary,
                              percentage_val = percentage[k], 
                              vec1 = basis_list[[k]]$rep1,
                              vec2 = basis_list[[k]]$rep2)
})

common_representation <- sapply(1:r, function(k){
  tiltedCCA:::.position_from_circle(circle_list[[k]], radian_vec[k])
})

common_score <- sapply(1:r, function(k){
  basis_list[[k]]$basis_mat %*% common_representation[,k]
})

common_mat <- tiltedCCA:::.convert_common_score_to_mat(common_score,
                                                       score_1,
                                                       score_2,
                                                       svd_1,
                                                       svd_2)

############################################

set.seed(10)
seurat_umap <- Seurat::RunUMAP(common_mat, 
                               assay = "SCT")
rownames(seurat_umap@cell.embeddings) <- rownames(score_1)

pbmc[["tmp"]] <- Seurat::CreateDimReducObject(seurat_umap@cell.embeddings,
                             assay = "SCT")

plot1 <- Seurat::DimPlot(pbmc, reduction = "tmp",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         raster=FALSE)
plot1 <- plot1 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nCommon subspace"))
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_tcca-umap_common_tmp.png"),
                plot1, device = "png", width = 6, height = 5, units = "in")

