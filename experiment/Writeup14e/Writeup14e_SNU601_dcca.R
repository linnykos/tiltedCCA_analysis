rm(list=ls())
load("../../../../out/Writeup14e/Writeup14e_SNU_preprocessed2.RData")
library(Seurat)
library(Signac)

mat_1 <- Matrix::t(SNU[["CNA"]]@counts)
mat_2 <- Matrix::t(SNU[["ATAC"]]@data)

zz <- svd(mat_1)$d
round(zz,2)

set.seed(10)
rank_1 <- 7; rank_2 <- 15; nn <- 30
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, 
                                      dims_1 = 1:rank_1, 
                                      dims_2 = 2:rank_2,
                                      discretization_gridsize = 17,
                                      metacell_clustering = NA, 
                                      num_neigh = nn, 
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      fix_tilt_perc = 0.5, 
                                      verbose = T) 

names(dcca_res)
dcca_res$df_percentage
dcca_res$cca_obj
dcca_res$tilt_perc
table(dcca_res$metacell_clustering, SNU@meta.data$clone)
head(svd(dcca_res$distinct_score_1)$d)

#############

print(paste0(Sys.time(),": Copy number frNN"))
set.seed(10)
cna_frnn <- multiomicCCA::construct_frnn(dcca_res, 
                                         nn = nn, 
                                         membership_vec = NA,
                                         data_1 = T, 
                                         data_2 = F,
                                         normalization_type = "none",
                                         radius_quantile = 0.5, 
                                         symmetrize = F, 
                                         bool_matrix = T, 
                                         verbose = T)

print(paste0(Sys.time(),": Copy number embedding"))
set.seed(10)
cna_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, 
                                                 nn = nn, 
                                                 data_1 = T, 
                                                 data_2 = F,
                                                 c_g = cna_frnn$c_g, 
                                                 d_g = cna_frnn$d_g, 
                                                 only_embedding = T, 
                                                 verbose = T, 
                                                 sampling_type = "adaptive_gaussian",
                                                 keep_nn = T)

print(paste0(Sys.time(),": ATAC frNN"))
set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res, 
                                         nn = nn, 
                                         membership_vec = NA,
                                         data_1 = F, 
                                         data_2 = T,
                                         normalization_type = "signac_itself",
                                         radius_quantile = 0.5,
                                         bool_matrix = T, symmetrize = F, verbose = T)

print(paste0(Sys.time(),": ATAC embedding"))
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
                                                 keep_nn = T)

print(paste0(Sys.time(),": Combining embeddings"))
set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, 
                                         g_1 = cna_frnn$c_g,
                                         g_2 = atac_frnn$c_g,
                                         nn = nn, 
                                         verbose = 1)
combined_g <- SeuratObject::as.Graph(combined_g)
set.seed(10)
combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings

save(SNU, dcca_res, cna_frnn, cna_embeddings,
     atac_frnn, atac_embeddings, combined_common_umap,
     file = paste0("../../../../out/Writeup14e/Writeup14e_SNU_dcca.RData"))
