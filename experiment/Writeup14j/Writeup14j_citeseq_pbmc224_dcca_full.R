rm(list=ls())
load("../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_preprocessed.RData")
library(Seurat)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()
n <- ncol(pbmc)

membership_vec <- as.factor(pbmc@meta.data$celltype.l2)
n <- length(membership_vec)
sort(table(membership_vec))
set.seed(10)
idx <- multiomicCCA::construct_celltype_subsample(membership_vec, 
                                                  min_subsample_cell = 5000)
keep_vec <- rep(0, n)
names(keep_vec) <- rownames(pbmc@meta.data)
keep_vec[idx] <- 1
pbmc$keep <- keep_vec
pbmc <- subset(pbmc, keep == 1)
pbmc[["rna.umap"]] <- NULL
pbmc[["adt.umap"]] <- NULL
pbmc[["wnn.umap"]] <- NULL

set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'pca', dims = 1:40, assay = 'SCT',
                        reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'apca', dims = 1:50, assay = 'ADT',
                        reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

############

n <- ncol(pbmc)
Seurat::DefaultAssay(pbmc) <- "SCT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:40)
set.seed(10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)
Seurat::DefaultAssay(pbmc) <- "ADT"
set.seed(10)
pbmc <- Seurat::FindNeighbors(pbmc, dims = 1:50, reduction = "apca")
set.seed(10)
pbmc <- Seurat::FindClusters(pbmc, resolution = 0.25)

###############

mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ADT"]]@scale.data)

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

rank_1 <- 40; rank_2 <- 50; nn <- 30
dims_1 <- 1:rank_1; dims_2 <- 1:rank_2

svd_1 <- multiomicCCA:::.svd_truncated(mat_1, K = rank_1, symmetric = F, rescale = F, 
                                       mean_vec = F, sd_vec = F, K_full_rank = F)
svd_2 <- multiomicCCA:::.svd_truncated(mat_2, K = rank_2, symmetric = F, rescale = F, 
                                       mean_vec = F, sd_vec = F, K_full_rank = F)

tmp <- multiomicCCA:::form_snns(num_neigh = 60,
                                svd_1 = svd_1,
                                svd_2 = svd_2,
                                bool_intersect = T,
                                distance_func = "cosine",
                                min_deg = 30,
                                verbose = T)
snn_mat_1 <- tmp$snn_mat_1; snn_mat_2 <- tmp$snn_mat_2
basis_1 <- multiomicCCA:::compute_laplacian_basis(snn_mat_1, k = 50, verbose = T)
basis_2 <- multiomicCCA:::compute_laplacian_basis(snn_mat_2, k = 50, verbose = T)


######################

clustering_1 <- factor(pbmc$SCT_snn_res.0.25)
clustering_2 <- factor(pbmc$ADT_snn_res.0.25)

# common neighbor
set.seed(10)
common_mat <- multiomicCCA:::common_neighborhood(snn_mat_1 = snn_mat_1, 
                                                 snn_mat_2 = snn_mat_2,
                                                 clustering_1 = clustering_1, 
                                                 clustering_2 = clustering_2,
                                                 num_neigh = 30,
                                                 verbose = T)
common_basis <- multiomicCCA:::compute_laplacian_basis(common_mat, k = 50, verbose = T)

save.image("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_full.RData")

set.seed(10)
dcca_res <- multiomicCCA::dcca_factor(mat_1b, mat_2b, 
                                      dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      fix_tilt_perc = F, 
                                      enforce_boundary = F,
                                      discretization_gridsize = 9,
                                      target_dimred = common_basis,
                                      verbose = T)
save.image("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_full.RData")

dcca_res2 <- multiomicCCA:::fine_tuning(dcca_res, 
                                        temp_path = "../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_full_tmp.RData",
                                        verbose = T)
dcca_res2$tilt_perc

save.image("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_full.RData")


