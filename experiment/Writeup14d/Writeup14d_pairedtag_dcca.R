rm(list=ls())
library(Seurat)
library(Signac)

histone_names <- c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3")
source("metacells.R")

for(i in 1:length(histone_names)){
  print(i)
  
  load(paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_",
              histone_names[i], ".RData"))
  
  # embedding_1 <- pairedtag[["pca"]]@cell.embeddings
  # embedding_2 <- pairedtag[["lsi"]]@cell.embeddings
  # 
  # set.seed(10)
  # meta_clustering <- compute_metacells(pairedtag,
  #                                      embedding_1 = embedding_1,
  #                                      embedding_2 = embedding_2)
  
  Seurat::DefaultAssay(pairedtag) <- "SCT"
  mat_1 <- Matrix::t(pairedtag[["SCT"]]@data)
  mat_2 <- Matrix::t(pairedtag[["DNA"]]@data)
  cell_name <- rownames(mat_1)
  membership_vec <- as.factor(pairedtag@meta.data$celltype)
  
  set.seed(10)
  rank_1 <- 50; rank_2 <- 25; nn <- 30
  dcca_res <- multiomicCCA::dcca_factor(mat_1, mat_2, 
                                        dims_1 = 1:rank_1, 
                                        dims_2 = 1:rank_2,
                                        meta_clustering = NA, 
                                        num_neigh = nn, 
                                        center_1 = T, center_2 = T,
                                        scale_1 = T, scale_2 = T,
                                        fix_distinct_perc = F, 
                                        verbose = T) 
  
  save(pairedtag, dcca_res, 
       file = paste0("../../../../out/Writeup14d/Writeup14d_pairedtag_",
                     histone_names[i], ".RData"))
  
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
  
  save(pairedtag, dcca_res, rna_frnn, rna_embeddings,
       file = paste0("../../../../out/Writeup14d/Writeup14d_pairedtag_",
                     histone_names[i], ".RData"))
  
  print(paste0(Sys.time(),": DNA frNN"))
  set.seed(10)
  dna_frnn <- multiomicCCA::construct_frnn(dcca_res, 
                                           nn = nn, 
                                           membership_vec = membership_vec,
                                           data_1 = F, 
                                           data_2 = T,
                                           normalization_type = "signac_itself",
                                           radius_quantile = 0.5,
                                           bool_matrix = T, symmetrize = F, verbose = T)
  
  print(paste0(Sys.time(),": DNA embedding"))
  set.seed(10)
  dna_embeddings <- multiomicCCA::plot_embeddings2(dcca_res, 
                                                   nn = nn, 
                                                   data_1 = F, 
                                                   data_2 = T,
                                                   c_g = dna_frnn$c_g, 
                                                   d_g = dna_frnn$d_g, 
                                                   only_embedding = T,
                                                   verbose = T, 
                                                   sampling_type = "adaptive_gaussian",
                                                   keep_nn = T)
  
  save(pairedtag, dcca_res, rna_frnn, rna_embeddings,
       dna_frnn, dna_embeddings,
       file = paste0("../../../../out/Writeup14d/Writeup14d_pairedtag_",
                     histone_names[i], ".RData"))
  
  print(paste0(Sys.time(),": Combining embeddings"))
  set.seed(10)
  combined_g <- multiomicCCA::combine_frnn(dcca_res, 
                                           g_1 = rna_frnn$c_g,
                                           g_2 = dna_frnn$c_g,
                                           nn = nn, 
                                           verbose = 2)
  combined_g <- SeuratObject::as.Graph(combined_g)
  set.seed(10)
  combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings
  
  save(pairedtag, dcca_res, rna_frnn, rna_embeddings,
       dna_frnn, dna_embeddings, combined_common_umap,
       file = paste0("../../../../out/Writeup14d/Writeup14d_pairedtag_",
                     histone_names[i], ".RData"))
  
}