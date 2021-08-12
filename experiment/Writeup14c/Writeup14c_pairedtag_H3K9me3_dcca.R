rm(list=ls())
load("../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_H3K9me3.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(pairedtag) <- "SCT"
mat_1 <- Matrix::t(pairedtag[["SCT"]]@data[Seurat::VariableFeatures(pairedtag),])

Seurat::DefaultAssay(pairedtag) <- "DNA"
mat_2 <- Matrix::t(pairedtag[["DNA"]]@data[Seurat::VariableFeatures(pairedtag),])

set.seed(10)
rank_1 <- 50; rank_2 <- 50; nn <- 30
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 1:rank_2,
                                      meta_clustering = NA, num_neigh = nn,
                                      center_1 = F, center_2 = F,
                                      scale_1 = F, scale_2 = F,
                                      fix_distinct_perc = F, 
                                      verbose = T) 

save(date_of_run, session_info, dcca_res,
     file = "../../../../out/Writeup14c/Writeup14c_pairedtag_H3K9me3_dcca.RData")

###############################

pairedtag2 <- Seurat::CreateSeuratObject(counts = Matrix::t(mat_1[,1:10]))
pairedtag2[["celltype"]] <- pairedtag@meta.data$celltype
membership_vec <- as.factor(pairedtag@meta.data$celltype)
rm(list = c("pairedtag", "mat_1", "mat_2")); gc()

set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                         data_1 = T, data_2 = F,
                                         radius_quantile = 0.5, symmetrize = F, 
                                         bool_matrix = T, verbose = T)

#compute all rna basis vectors
k_max <- 200
c_eig <- multiomicCCA::compute_laplacian(rna_frnn$c_g, k_max = k_max, rowname_vec = colnames(pairedtag2), 
                                         colname_vec = paste0("clap_", 1:k_max), verbose = F)
d_eig <- multiomicCCA::compute_laplacian(rna_frnn$d_g, k_max = k_max, rowname_vec = colnames(pairedtag2), 
                                         colname_vec = paste0("dlap_", 1:k_max), verbose = F)
e_eig <- multiomicCCA::compute_laplacian(rna_frnn$e_g, k_max = k_max, rowname_vec = colnames(pairedtag2), 
                                         colname_vec = paste0("elap_", 1:k_max), verbose = F)

pairedtag2[["clap"]] <- Seurat::CreateDimReducObject(embedding = c_eig, key = "clap", assay = "RNA")
pairedtag2[["dlap"]] <- Seurat::CreateDimReducObject(embedding = d_eig, key = "dlap", assay = "RNA")
pairedtag2[["elap"]] <- Seurat::CreateDimReducObject(embedding = e_eig, key = "elap", assay = "RNA")

set.seed(10)
rna_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = T, data_2 = F,
                                                 c_g = rna_frnn$c_g, d_g = rna_frnn$d_g, 
                                                 only_embedding = T, verbose = T, 
                                                 sampling_type = "adaptive_gaussian",
                                                 metric = "euclidean",
                                                 keep_nn = T)
pairedtag2[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "UMAP", assay = "RNA")
pairedtag2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "UMAP", assay = "RNA")
pairedtag2[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "UMAP", assay = "RNA")

save(date_of_run, session_info, dcca_res, pairedtag2, rna_frnn,
     file = "../../../../out/Writeup14c/Writeup14c_pairedtag_H3K9me3_dcca.RData")
print("Finished RNA embedings")

###########################

set.seed(10)
histone_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                             data_1 = F, data_2 = T,
                                             radius_quantile = 0.5,
                                             bool_matrix = T, symmetrize = F, verbose = T)

#compute all the degree vectors
k_max <- 200
c_eig2 <- multiomicCCA::compute_laplacian(histone_frnn$c_g, k_max = k_max, rowname_vec = colnames(pairedtag2), 
                                          colname_vec = paste0("clap_", 1:k_max), verbose = F)
d_eig2 <- multiomicCCA::compute_laplacian(histone_frnn$d_g, k_max = k_max, rowname_vec = colnames(pairedtag2), 
                                          colname_vec = paste0("dlap_", 1:k_max), verbose = F)
e_eig2 <- multiomicCCA::compute_laplacian(histone_frnn$e_g, k_max = k_max, rowname_vec = colnames(pairedtag2), 
                                          colname_vec = paste0("elap_", 1:k_max), verbose = F)

pairedtag2[["clap2"]] <- Seurat::CreateDimReducObject(embedding = c_eig2, key = "clap", assay = "RNA")
pairedtag2[["dlap2"]] <- Seurat::CreateDimReducObject(embedding = d_eig2, key = "dlap", assay = "RNA")
pairedtag2[["elap2"]] <- Seurat::CreateDimReducObject(embedding = e_eig2, key = "elap", assay = "RNA")

set.seed(10)
histone_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = F, data_2 = T,
                                                     c_g = histone_frnn$c_g, d_g = histone_frnn$d_g, 
                                                     only_embedding = T, verbose = T, 
                                                     sampling_type = "adaptive_gaussian",
                                                     metric = "euclidean",
                                                     keep_nn = T)
pairedtag2[["common2"]] <- Seurat::CreateDimReducObject(embedding = histone_embeddings[[1]], key = "UMAP", assay = "RNA")
pairedtag2[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = histone_embeddings[[2]], key = "UMAP", assay = "RNA")
pairedtag2[["everything2"]] <- Seurat::CreateDimReducObject(embedding = histone_embeddings[[3]], key = "UMAP", assay = "RNA")

save(date_of_run, session_info, dcca_res, pairedtag2, rna_frnn, histone_frnn,
     file = "../../../../out/Writeup14c/Writeup14c_pairedtag_H3K9me3_dcca.RData")
print("Finished histone embedings")

##########################

set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, g_1 = rna_frnn$c_g,
                                         g_2 = histone_frnn$c_g, nn = nn, 
                                         verbose = 2)
combined_g <- SeuratObject::as.Graph(combined_g)
set.seed(10)
combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings
pairedtag2[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, 
                                                         key = "UMAP", assay = "RNA")

save(date_of_run, session_info, dcca_res, pairedtag2, rna_frnn, histone_frnn,
     file = "../../../../out/Writeup14c/Writeup14c_pairedtag_H3K9me3_dcca.RData")

set.seed(10)
both_embeddings <- multiomicCCA::plot_embeddings(dcca_res, membership_vec = membership_vec,
                                                 data_1 = T, data_2 = T, add_noise = F,
                                                 only_embedding = T)
pairedtag2[["both"]] <- Seurat::CreateDimReducObject(embedding = both_embeddings$everything, 
                                                     key = "UMAP", assay = "RNA")

save(date_of_run, session_info, dcca_res, pairedtag2, rna_frnn, histone_frnn,
     file = "../../../../out/Writeup14c/Writeup14c_pairedtag_H3K9me3_dcca.RData")
print("Finished joint embedings")

###########################

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
mat_2_denoised <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2
rm(list = c("dcca_decomp")); gc()

p1 <- ncol(mat_1_denoised)
gene_smoothed <- lapply(1:p1, function(j){
  if(j %% floor(p1/10) == 0) cat('*')
  
  c_res <- compute_smooth_signal(mat_1_denoised[,j], c_eig)
  d_res <- compute_smooth_signal(mat_1_denoised[,j], d_eig)
  e_res <- compute_smooth_signal(mat_1_denoised[,j], e_eig)
  
  list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
       d_variance = d_res$variance, d_r2 = d_res$r_squared,
       e_variance = e_res$variance, e_r2 = e_res$r_squared)
})

save(date_of_run, session_info, dcca_res, pairedtag2, rna_frnn, histone_frnn,
     gene_smoothed, file = "../../../../out/Writeup14c/Writeup14c_pairedtag_H3K9me3_dcca.RData")
print("Finished RNA variables")

###############################3

p2 <- ncol(mat_2_denoised)
histone_smoothed <- lapply(1:p2, function(j){
  if(j %% 10 == 0) print(j)
  
  c_res <- compute_smooth_signal(mat_2_denoised[,j], c_eig)
  d_res <- compute_smooth_signal(mat_2_denoised[,j], d_eig)
  e_res <- compute_smooth_signal(mat_2_denoised[,j], e_eig)
  
  list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
       d_variance = d_res$variance, d_r2 = d_res$r_squared,
       e_variance = e_res$variance, e_r2 = e_res$r_squared)
})

save(date_of_run, session_info, dcca_res, pairedtag2, rna_frnn, histone_frnn,
     gene_smoothed, histone_smoothed,
     file = "../../../../out/Writeup14c/Writeup14c_pairedtag_H3K9me3_dcca.RData")
print("Finished ATAC variables")


