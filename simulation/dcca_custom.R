dcca_custom <- function(dcca_res,
                        nn = 30){
  frnn_1 <- multiomicCCA::construct_frnn(dcca_res, 
                                         nn = nn, 
                                         membership_vec = NA,
                                         data_1 = T, 
                                         data_2 = F,
                                         normalization_type = "cosine_itself",
                                         radius_quantile = 0.5, 
                                         symmetrize = F, 
                                         bool_matrix = T, 
                                         verbose = F)
  
  dcca_embedding_1 <- multiomicCCA::plot_embeddings2(dcca_res, 
                                                     nn = nn, 
                                                     data_1 = T, 
                                                     data_2 = F,
                                                     c_g = frnn_1$c_g, 
                                                     d_g = frnn_1$d_g, 
                                                     only_embedding = T, 
                                                     verbose = F, 
                                                     sampling_type = "adaptive_gaussian",
                                                     keep_nn = T)
  
  frnn_2 <- multiomicCCA::construct_frnn(dcca_res, 
                                         nn = nn, 
                                         membership_vec = NA,
                                         data_1 = F, 
                                         data_2 = T,
                                         normalization_type = "cosine_itself",
                                         radius_quantile = 0.5,
                                         bool_matrix = T, 
                                         symmetrize = F, 
                                         verbose = F)
  
  
  dcca_embedding_2 <- multiomicCCA::plot_embeddings2(dcca_res, 
                                                     nn = nn, 
                                                     data_1 = F, 
                                                     data_2 = T,
                                                     c_g = frnn_2$c_g, 
                                                     d_g = frnn_2$d_g, 
                                                     only_embedding = T,
                                                     verbose = F, 
                                                     sampling_type = "adaptive_gaussian",
                                                     keep_nn = T)
  
  combined_g <- multiomicCCA::combine_frnn(dcca_res, 
                                           g_1 = frnn_1$c_g,
                                           g_2 = frnn_2$c_g,
                                           nn = nn, 
                                           verbose = 0)
  combined_g <- SeuratObject::as.Graph(combined_g)
  set.seed(10)
  combined_common_umap <- Seurat::RunUMAP(combined_g, assay = "RNA")@cell.embeddings
  
  list(common = combined_common_umap,
       distinct_1 = dcca_embedding_1[["distinct"]],
       distinct_2 = dcca_embedding_2[["distinct"]])
}
