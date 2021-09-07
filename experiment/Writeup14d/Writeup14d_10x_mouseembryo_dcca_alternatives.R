rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(multiomicCCA)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

source("frnn_alt.R")
source("plotting.R")

Seurat::DefaultAssay(mbrain) <- "SCT"
mat_1 <- Matrix::t(mbrain[["SCT"]]@data[Seurat::VariableFeatures(object = mbrain),])
Seurat::DefaultAssay(mbrain) <- "ATAC"
mat_2 <- Matrix::t(mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),])

head(rownames(mat_1)); head(colnames(mat_1))
head(rownames(mat_2)); head(colnames(mat_2))
dim(mat_1); dim(mat_2)
metadata <- mbrain@meta.data

cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Radial glia", 
                                                 "Forebrain GABAergic", "Neuroblast", 
                                                 "Glioblast", "Cortical or hippocampal glutamatergic"))
####################

# modification 1: removing both mean and standard deviation
mat_1b <- mat_1[cell_idx,]
mean_val_1 <- mean(mat_1b@x)
row_vec_1 <- sparseMatrixStats::rowSums2(mat_1b)
tmp <- Matrix::Diagonal(nrow(mat_1b), 1/row_vec_1)
mat_1b <- tmp %*% mat_1b 
mat_1b <- mat_1b * mean_val_1 / mean(mat_1b@x)

mat_2b <- mat_2[cell_idx,]
mean_val_2 <- mean(mat_2b@x)
row_vec_2 <- sparseMatrixStats::rowSums2(mat_2b)
tmp <- Matrix::Diagonal(nrow(mat_2b), 1/row_vec_2)
mat_2b <- tmp %*% mat_2b 
mat_2b <- mat_2b * mean_val_2 / mean(mat_2b@x)

set.seed(10)
rank_1 <- 30; rank_2 <- 50; nn <- 15
dcca_res <- multiomicCCA::dcca_factor(mat_1b, mat_2b, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = T,
                                      meta_clustering = NA, num_neigh = nn, 
                                      fix_distinct_perc = F, verbose = T) 
print("Finished Tilted-CCA")

mbrain2 <- Seurat::CreateSeuratObject(counts = Matrix::t(mat_1b[,1:10]))
mbrain2[["celltype"]] <- metadata$label_Savercat[cell_idx]
membership_vec <- as.factor(metadata$label_Savercat[cell_idx])
rm(list = c("mbrain", "mat_1", "mat_2")); gc()

set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                         data_1 = T, data_2 = F,
                                         radius_quantile = 0.5, symmetrize = F, 
                                         bool_matrix = T, verbose = T)

set.seed(10)
rna_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = T, data_2 = F,
                                                 c_g = rna_frnn$c_g, d_g = rna_frnn$d_g, 
                                                 only_embedding = T, verbose = T, 
                                                 sampling_type = "adaptive_gaussian",
                                                 keep_nn = T)
mbrain2[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "UMAP", assay = "RNA")
mbrain2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "UMAP", assay = "RNA")
mbrain2[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "UMAP", assay = "RNA")

set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                          data_1 = F, data_2 = T,
                                          radius_quantile = 0.5,
                                          bool_matrix = T, symmetrize = F, verbose = T)

set.seed(10)
atac_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = F, data_2 = T,
                                                  c_g = atac_frnn$c_g, d_g = atac_frnn$d_g, 
                                                  only_embedding = T, verbose = T, 
                                                  sampling_type = "adaptive_gaussian",
                                                  keep_nn = T)
mbrain2[["common2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[1]], key = "UMAP", assay = "RNA")
mbrain2[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[2]], key = "UMAP", assay = "RNA")
mbrain2[["everything2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[3]], key = "UMAP", assay = "RNA")

set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, g_1 = rna_frnn$c_g,
                                         g_2 = atac_frnn$c_g, nn = nn, 
                                         verbose = 1)
combined_g <- SeuratObject::as.Graph(combined_g)
set.seed(10)
combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings
mbrain2[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, 
                                                      key = "UMAP", assay = "RNA")

set.seed(10)
both_embeddings <- multiomicCCA::plot_embeddings(dcca_res, membership_vec = membership_vec,
                                                 data_1 = T, data_2 = T, add_noise = F,
                                                 only_embedding = T)
mbrain2[["both"]] <- Seurat::CreateDimReducObject(embedding = both_embeddings$everything, 
                                                  key = "UMAP", assay = "RNA")

title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

make_all_embedding_plots(mbrain2, 
                         title_all = "Mouse Embryo (Row-rescale)",
                         title_vec = title_vec, 
                         main_vec = main_vec, 
                         file_prefix = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_dcca",
                         file_suffix = "meansdrowrescale")

###############################

# modification 1: removing both mean and standard deviation
mat_1b <- mat_1[cell_idx,]
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2[cell_idx,]
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

set.seed(10)
rank_1 <- 30; rank_2 <- 50; nn <- 15
dcca_res <- multiomicCCA::dcca_factor(mat_1b, mat_2b, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      center_1 = T, center_2 = T,
                                      meta_clustering = NA, num_neigh = nn, 
                                      fix_distinct_perc = F, verbose = T) 
print("Finished Tilted-CCA")

mbrain2 <- Seurat::CreateSeuratObject(counts = Matrix::t(mat_1b[,1:10]))
mbrain2[["celltype"]] <- metadata$label_Savercat[cell_idx]
membership_vec <- as.factor(metadata$label_Savercat[cell_idx])
rm(list = c("mbrain", "mat_1", "mat_2")); gc()

set.seed(10)
rna_frnn <- construct_frnn_alt(dcca_res, 
                               nn = nn, 
                               membership_vec = membership_vec,
                               data_1 = T, 
                               data_2 = F,
                               radius_quantile = 0.5, 
                               symmetrize = F, 
                               bool_matrix = T, 
                               verbose = T,
                               normalization_type = "everything")


set.seed(10)
rna_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, 
                                                 nn = nn, 
                                                 data_1 = T, 
                                                 data_2 = F,
                                                 c_g = rna_frnn$c_g, 
                                                 d_g = rna_frnn$d_g, 
                                                 only_embedding = T, 
                                                 verbose = T, 
                                                 sampling_type = "adaptive_gaussian",
                                                 keep_nn = T,
                                                 metric = "cosine")
mbrain2[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "common_", assay = "RNA")
mbrain2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "distinct_", assay = "RNA")
mbrain2[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "everything_", assay = "RNA")

set.seed(10)
atac_frnn <- construct_frnn_alt(dcca_res, 
                                nn = nn, 
                                membership_vec = membership_vec,
                                data_1 = F, 
                                data_2 = T,
                                radius_quantile = 0.5,
                                bool_matrix = T, 
                                symmetrize = F, 
                                verbose = T,
                                normalization_type = "everything")

set.seed(10)
atac_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, 
                                                  nn = nn, 
                                                  data_1 = F, 
                                                  data_2 = T,
                                                  c_g = atac_frnn$c_g,
                                                  d_g = atac_frnn$d_g, 
                                                  only_embedding = T, 
                                                  verbose = T, 
                                                  sampling_type = "adaptive_gaussian",
                                                  keep_nn = T,
                                                  metric = "cosine")
mbrain2[["common2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[1]], key = "common2_", assay = "RNA")
mbrain2[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[2]], key = "distinct2_", assay = "RNA")
mbrain2[["everything2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[3]], key = "everything2_", assay = "RNA")

set.seed(10)
combined_g <- combine_frnn_alt(dcca_res, 
                               g_1 = rna_frnn$c_g,
                               g_2 = atac_frnn$c_g, 
                               nn = nn, 
                               verbose = 1)
combined_g <- SeuratObject::as.Graph(combined_g)
set.seed(10)
combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings
mbrain2[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, 
                                                      key = "combined_", assay = "RNA")

set.seed(10)
both_embeddings <- multiomicCCA::plot_embeddings(dcca_res, membership_vec = membership_vec,
                                                 data_1 = T, data_2 = T, add_noise = F,
                                                 only_embedding = T)
mbrain2[["both"]] <- Seurat::CreateDimReducObject(embedding = both_embeddings$everything, 
                                                  key = "both_", assay = "RNA")

title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

make_all_embedding_plots(mbrain2, 
                         title_all = "Mouse Embryo (L2norm-everything)",
                         title_vec = title_vec, 
                         main_vec = main_vec, 
                         file_prefix = "../../../../out/figures/Writeup14d/Writeup14d_10x_mouseembryo_dcca",
                         file_suffix = "_l2everything")

