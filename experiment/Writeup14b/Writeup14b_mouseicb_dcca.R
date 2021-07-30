rm(list=ls())
load("../../../../out/Writeup14/Writeup14_mouseicb_preprocess.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

membership_vec <- as.character(myeloid@meta.data$sample)
membership_vec <- sapply(membership_vec, function(x){
        paste0(strsplit(x, split = "_")[[1]][1:2], collapse = "_")
})
names(membership_vec) <- NULL
membership_vec <- as.factor(membership_vec)
table(membership_vec)
range(table(membership_vec))

mat_1 <- t(myeloid[["SCT"]]@scale.data)

Seurat::DefaultAssay(myeloid) <- "ATAC"
myeloid <- Seurat::ScaleData(myeloid)
mat_2 <- t(myeloid[["ATAC"]]@scale.data)

set.seed(10)
rank_1 <- 50; rank_2 <- 50; nn <- 30
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, dims_1 = 1:rank_1, dims_2 = 2:rank_2,
                                      meta_clustering = NA, num_neigh = nn,
                                      fix_distinct_perc = F, verbose = T) 

rm(list = c("myeloid", "mat_1", "mat_2")); gc()

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
mat_2_denoised <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2

save(date_of_run, session_info, dcca_res,
     file = "../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

rm(list = "dcca_decomp"); gc()

###############

myeloid2 <- Seurat::CreateSeuratObject(counts = t(mat_1_denoised))
myeloid2[["celltype"]] <- membership_vec

title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                         data_1 = T, data_2 = F,
                                         radius_quantile = 0.5, symmetrize = F, 
                                         bool_matrix = T, verbose = T)

#compute all rna basis vectors
k_max <- 200
c_eig <- multiomicCCA::compute_laplacian(rna_frnn$c_g, k_max = k_max, rowname_vec = colnames(myeloid2), 
                                         colname_vec = paste0("clap_", 1:k_max), verbose = F)
d_eig <- multiomicCCA::compute_laplacian(rna_frnn$d_g, k_max = k_max, rowname_vec = colnames(myeloid2), 
                                         colname_vec = paste0("dlap_", 1:k_max), verbose = F)
e_eig <- multiomicCCA::compute_laplacian(rna_frnn$e_g, k_max = k_max, rowname_vec = colnames(myeloid2), 
                                         colname_vec = paste0("elap_", 1:k_max), verbose = F)

myeloid2[["clap"]] <- Seurat::CreateDimReducObject(embedding = c_eig, key = "clap", assay = "RNA")
myeloid2[["dlap"]] <- Seurat::CreateDimReducObject(embedding = d_eig, key = "dlap", assay = "RNA")
myeloid2[["elap"]] <- Seurat::CreateDimReducObject(embedding = e_eig, key = "elap", assay = "RNA")

set.seed(10)
rna_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = T, data_2 = F,
                                                 c_g = rna_frnn$c_g, d_g = rna_frnn$d_g, 
                                                 only_embedding = T, verbose = T, 
                                                 sampling_type = "adaptive_gaussian",
                                                 keep_nn = T)
myeloid2[["common"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[1]], key = "UMAP", assay = "RNA")
myeloid2[["distinct"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[2]], key = "UMAP", assay = "RNA")
myeloid2[["everything"]] <- Seurat::CreateDimReducObject(embedding = rna_embeddings[[3]], key = "UMAP", assay = "RNA")

save(date_of_run, session_info, dcca_res, myeloid2, rna_frnn,
     file = "../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

################

set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res, nn = nn, membership_vec = membership_vec,
                                             data_1 = F, data_2 = T,
                                             radius_quantile = 0.5,
                                             bool_matrix = T, symmetrize = F, verbose = T)

#compute all the degree vectors
k_max <- 200
c_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$c_g, k_max = k_max, rowname_vec = colnames(myeloid2), 
                                          colname_vec = paste0("clap_", 1:k_max), verbose = F)
d_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$d_g, k_max = k_max, rowname_vec = colnames(myeloid2), 
                                          colname_vec = paste0("dlap_", 1:k_max), verbose = F)
e_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$e_g, k_max = k_max, rowname_vec = colnames(myeloid2), 
                                          colname_vec = paste0("elap_", 1:k_max), verbose = F)

myeloid2[["clap2"]] <- Seurat::CreateDimReducObject(embedding = c_eig2, key = "clap", assay = "RNA")
myeloid2[["dlap2"]] <- Seurat::CreateDimReducObject(embedding = d_eig2, key = "dlap", assay = "RNA")
myeloid2[["elap2"]] <- Seurat::CreateDimReducObject(embedding = e_eig2, key = "elap", assay = "RNA")

set.seed(10)
atac_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, nn = nn, data_1 = F, data_2 = T,
                                                     c_g = atac_frnn$c_g, d_g = atac_frnn$d_g, 
                                                     only_embedding = T, verbose = T, 
                                                     sampling_type = "adaptive_gaussian",
                                                     keep_nn = T)
myeloid2[["common2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[1]], key = "UMAP", assay = "RNA")
myeloid2[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[2]], key = "UMAP", assay = "RNA")
myeloid2[["everything2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[3]], key = "UMAP", assay = "RNA")

save(date_of_run, session_info, dcca_res, myeloid2, rna_frnn, atac_frnn,
     file = "../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

##########################

set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, g_1 = rna_frnn$c_g,
                                         g_2 = atac_frnn$c_g, nn = nn, 
                                         verbose = 2)
combined_g <- SeuratObject::as.Graph(combined_g)
set.seed(10)
combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings
myeloid2[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, key = "UMAP", assay = "RNA")

save(date_of_run, session_info, dcca_res, myeloid2, rna_frnn, atac_frnn, combined_g,
     file = "../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

set.seed(10)
both_embeddings <- multiomicCCA::plot_embeddings(dcca_res, membership_vec = membership_vec,
                                                 data_1 = T, data_2 = T, add_noise = F,
                                                 only_embedding = T)
myeloid2[["both"]] <- Seurat::CreateDimReducObject(embedding = both_embeddings$everything, key = "UMAP", assay = "RNA")

save(date_of_run, session_info, dcca_res, myeloid2, rna_frnn, atac_frnn, combined_g,
     file = "../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

#################

print("Starting RNA smooth")
p1 <- ncol(mat_1_denoised)
gene_smoothed <- lapply(1:p1, function(j){
        print(j)
        
        c_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,j], c_eig)
        d_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,j], d_eig)
        e_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,j], e_eig)
        
        list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
             d_variance = d_res$variance, d_r2 = d_res$r_squared,
             e_variance = e_res$variance, e_r2 = e_res$r_squared)
})

save(date_of_run, session_info, dcca_res, myeloid2, rna_frnn, atac_frnn, combined_g,
     gene_smoothed, 
     file = "../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

print("Starting ATAC smooth")
p2 <- ncol(mat_2_denoised)
atac_smoothed <- lapply(1:p2, function(j){
        print(j)
        
        c_res <- multiomicCCA::compute_smooth_signal(mat_2_denoised[,j], c_eig2)
        d_res <- multiomicCCA::compute_smooth_signal(mat_2_denoised[,j], d_eig2)
        e_res <- multiomicCCA::compute_smooth_signal(mat_2_denoised[,j], e_eig2)
        
        list(c_variance = c_res$variance, c_r2 = c_res$r_squared,
             d_variance = d_res$variance, d_r2 = d_res$r_squared,
             e_variance = e_res$variance, e_r2 = e_res$r_squared)
})

save(date_of_run, session_info, dcca_res, myeloid2, rna_frnn, atac_frnn, combined_g,
     gene_smoothed, atac_smoothed,
     file = "../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")
