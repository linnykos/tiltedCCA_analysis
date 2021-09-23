rm(list=ls())
load("../../../../out/Writeup11/Writeup11_citeseq_pbmc228_preprocessed2.RData")

library(Seurat)
library(multiomicCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

mat_1 <- t(pbmc[["SCT"]]@scale.data)
mat_2 <- t(pbmc[["ADT"]]@scale.data)
dim(mat_1); dim(mat_2)

set.seed(10)
idx <- multiomicCCA::construct_celltype_subsample(membership_vec, min_subsample_cell = 3000)
mat_1 <- mat_1[idx,]; mat_2 <- mat_2[idx,]
cell_name <- rownames(mat_1)
membership_vec <- as.factor(pbmc@meta.data$celltype.l2[idx])

print(paste0(Sys.time(),": Applying D-CCA"))
set.seed(10)
rank_1 <- 40; rank_2 <- 50; nn <- 30
dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, 
                                      dims_1 = 1:rank_1, 
                                      dims_2 = 1:rank_2,
                                      meta_clustering = NA, 
                                      num_neigh = nn, 
                                      cell_max = 15000,
                                      fix_distinct_perc = F, 
                                      verbose = T) 

rm(list = c("mat_1", "mat_2")); gc(T)

save.image("../../../../out/Writeup14d/Writeup14d_citeseq_pbmc228_dcca.RData")

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

save.image("../../../../out/Writeup14d/Writeup14d_citeseq_pbmc228_dcca.RData")

################

print(paste0(Sys.time(),": Protein frNN"))
set.seed(10)
protein_frnn <- multiomicCCA::construct_frnn(dcca_res, 
                                             nn = nn, 
                                             membership_vec = membership_vec,
                                             data_1 = F, 
                                             data_2 = T,
                                             normalization_type = "cosine_itself",
                                             radius_quantile = 0.5,
                                             bool_matrix = T, symmetrize = F, verbose = T)

#compute all the degree vectors
print(paste0(Sys.time(),": Protein eigenbasis"))
k_max <- 50
c_eig2 <- multiomicCCA::compute_laplacian(protein_frnn$c_g, 
                                          k_max = k_max, 
                                          rowname_vec = cell_name, 
                                          colname_vec = paste0("clap_", 1:k_max), verbose = F)
d_eig2 <- multiomicCCA::compute_laplacian(protein_frnn$d_g,
                                          k_max = k_max, 
                                          rowname_vec = cell_name, 
                                          colname_vec = paste0("dlap_", 1:k_max), verbose = F)
e_eig2 <- multiomicCCA::compute_laplacian(protein_frnn$e_g, 
                                          k_max = k_max, 
                                          rowname_vec = cell_name, 
                                          colname_vec = paste0("elap_", 1:k_max), verbose = F)

print(paste0(Sys.time(),": Protein embedding"))
set.seed(10)
protein_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, 
                                                     nn = nn, 
                                                     data_1 = F, 
                                                     data_2 = T,
                                                     c_g = protein_frnn$c_g, 
                                                     d_g = protein_frnn$d_g, 
                                                     only_embedding = T,
                                                     verbose = T, 
                                                     sampling_type = "adaptive_gaussian",
                                                     keep_nn = T)

save.image("../../../../out/Writeup14d/Writeup14d_citeseq_pbmc228_dcca.RData")

##########################

print(paste0(Sys.time(),": Combining embeddings"))
set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, 
                                         g_1 = rna_frnn$c_g,
                                         g_2 = protein_frnn$c_g,
                                         nn = nn, 
                                         verbose = 2)
combined_g <- SeuratObject::as.Graph(combined_g)
set.seed(10)
combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings

save.image("../../../../out/Writeup14d/Writeup14d_citeseq_pbmc228_dcca.RData")

print("Finished")