rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(EnsDb.Hsapiens.v86); 
library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
library(dplyr); library(ggplot2); library(multiomicCCA)

load("../../../../out/Writeup10/Writeup10_10x_pbmc_preprocess4.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

source("metacells.R")

tab <- table(pbmc@meta.data$predicted.id)
rm_idx <- which(pbmc@meta.data$predicted.id %in% names(tab)[tab <= 5])
keep_vec <- rep(1, ncol(pbmc))
keep_vec[rm_idx] <- 0
pbmc[["keep"]] <- keep_vec
pbmc <- subset(pbmc, keep == 1)

embedding_1 <- pbmc[["pca"]]@cell.embeddings
embedding_2 <- pbmc[["lsi"]]@cell.embeddings

set.seed(10)
meta_clustering <- compute_metacells(pbmc,
                                     embedding_1 = embedding_1,
                                     embedding_2 = embedding_2)

mat_1 <- Matrix::t(pbmc[["SCT"]]@data)
mat_2 <- Matrix::t(pbmc[["ATAC"]]@data)
cell_name <- rownames(mat_1)
membership_vec <- as.factor(pbmc@meta.data$predicted.id)

set.seed(10)
rank_1 <- 40; rank_2 <- 50; nn <- 30
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, 
                                      dims_1 = 1:rank_1, 
                                      dims_2 = 2:rank_2,
                                      metacell_clustering = meta_clustering, 
                                      num_neigh = nn, 
                                      center_1 = T, center_2 = F,
                                      scale_1 = T, scale_2 = F,
                                      fix_tilt_perc = F, 
                                      verbose = T) 


rm(list = c("mat_1", "mat_2")); gc(T)

save.image("../../../../out/Writeup14d/Writeup14d_10x_pbmc_dcca.RData")

##################

print(paste0(Sys.time(),": RNA frNN"))
set.seed(10)
rna_frnn <- multiomicCCA::construct_frnn(dcca_res, 
                                         nn = nn, 
                                         membership_vec = membership_vec,
                                         data_1 = T, 
                                         data_2 = F,
                                         normalization_type = "cosine_itself",
                                         radius_quantile = 0.5, 
                                         symmetrize = F, 
                                         bool_matrix = T, 
                                         verbose = T)

#compute all rna basis vectors
print(paste0(Sys.time(),": RNA eigenbasis"))
k_max <- 50
c_eig <- multiomicCCA::compute_laplacian(rna_frnn$c_g, 
                                         k_max = k_max, 
                                         rowname_vec = cell_name, 
                                         colname_vec = paste0("clap_", 1:k_max), verbose = F)
d_eig <- multiomicCCA::compute_laplacian(rna_frnn$d_g, 
                                         k_max = k_max,
                                         rowname_vec = cell_name, 
                                         colname_vec = paste0("dlap_", 1:k_max), verbose = F)
e_eig <- multiomicCCA::compute_laplacian(rna_frnn$e_g, 
                                         k_max = k_max, 
                                         rowname_vec = cell_name, 
                                         colname_vec = paste0("elap_", 1:k_max), verbose = F)

print(paste0(Sys.time(),": RNA embedding"))
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
                                                 keep_nn = T)

save.image("../../../../out/Writeup14d/Writeup14d_10x_pbmc_dcca.RData")

################

print(paste0(Sys.time(),": ATAC frNN"))
set.seed(10)
atac_frnn <- multiomicCCA::construct_frnn(dcca_res, 
                                          nn = nn, 
                                          membership_vec = membership_vec,
                                          data_1 = F, 
                                          data_2 = T,
                                          normalization_type = "signac_itself",
                                          radius_quantile = 0.5,
                                          bool_matrix = T, symmetrize = F, verbose = T)

#compute all the degree vectors
print(paste0(Sys.time(),": ATAC eigenbasis"))
k_max <- 50
c_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$c_g, 
                                          k_max = k_max, 
                                          rowname_vec = cell_name, 
                                          colname_vec = paste0("clap_", 1:k_max), verbose = F)
d_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$d_g,
                                          k_max = k_max, 
                                          rowname_vec = cell_name, 
                                          colname_vec = paste0("dlap_", 1:k_max), verbose = F)
e_eig2 <- multiomicCCA::compute_laplacian(atac_frnn$e_g, 
                                          k_max = k_max, 
                                          rowname_vec = cell_name, 
                                          colname_vec = paste0("elap_", 1:k_max), verbose = F)

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

save.image("../../../../out/Writeup14d/Writeup14d_10x_pbmc_dcca.RData")

##########################

print(paste0(Sys.time(),": Combining embeddings"))
set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, 
                                         g_1 = rna_frnn$c_g,
                                         g_2 = atac_frnn$c_g,
                                         nn = nn, 
                                         verbose = 2)
combined_g <- SeuratObject::as.Graph(combined_g)
set.seed(10)
combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings

save.image("../../../../out/Writeup14d/Writeup14d_10x_pbmc_dcca.RData")

print("Finished")

