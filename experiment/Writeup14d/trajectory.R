compute_principal_ordering <- function(obj,
                                       clustering, 
                                       selected_clusters,
                                       selected_clusters_ordering,
                                       target_dimension = 5){
  
  # step 1: extract the embeddings
  tmp <- multiomicCCA:::.prepare_embeddings(obj, 
                                            data_1 = T, 
                                            data_2 = F, 
                                            add_noise = F, 
                                            center = F, 
                                            renormalize = F)
  rna_embedding <- tmp[["everything"]]
  l2_vec <- apply(rna_embedding, 1, multiomicCCA:::.l2norm)
  rna_embedding <- multiomicCCA:::.mult_vec_mat(1/l2_vec, rna_embedding)
  rna_embedding <- rna_embedding[clustering %in% selected_clusters,]
  rna_embedding_meta <-.form_metacell_matrix(rna_embedding, clustering[clustering %in% selected_clusters],
                                             func = mean)

  tmp <- multiomicCCA:::.prepare_embeddings(obj, 
                                            data_1 = F, 
                                            data_2 = T, 
                                            add_noise = F, 
                                            center = F, 
                                            renormalize = F)
  atac_embedding <- tmp[["everything"]]
  center_vec <- matrixStats::colMeans2(atac_embedding)
  sd_vec <- matrixStats::colSds(atac_embedding)
  for(j in 1:ncol(atac_embedding_meta)){
    atac_embedding[,j] <- (atac_embedding[,j]-center_vec[j])/sd_vec[j]
  }
  l2_vec <- apply(atac_embedding, 1, multiomicCCA:::.l2norm)
  atac_embedding <- multiomicCCA:::.mult_vec_mat(1/l2_vec, atac_embedding)
  atac_embedding <- atac_embedding[clustering %in% selected_clusters,]
  atac_embedding_meta <-.form_metacell_matrix(atac_embedding, clustering[clustering %in% selected_clusters],
                                              func = mean)
  
  # # step 2: reduce the dimensionality of both dimensions
  svd_res <- svd(rna_embedding_meta)
  rna_embedding_meta <- rna_embedding_meta %*% svd_res$v[,1:target_dimension,drop = F]
  rna_embedding <- rna_embedding %*% svd_res$v[,1:target_dimension,drop = F]

  svd_res <- svd(atac_embedding_meta)
  atac_embedding_meta <- atac_embedding_meta %*% svd_res$v[,1:target_dimension,drop = F]
  atac_embedding <- atac_embedding %*% svd_res$v[,1:target_dimension,drop = F]

  mat_meta <- 100*cbind(rna_embedding_meta, atac_embedding_meta)
  mat <- 100*cbind(rna_embedding, atac_embedding)
  
  # step 3a: reorder mat_meta
  mat_meta <- mat_meta[sapply(selected_clusters, function(x){
    which(rownames(mat_meta) == x)
  }),]
  mat_meta <- mat_meta[sapply(sort(unique(selected_clusters_ordering)), function(x){
    which(selected_clusters_ordering == x)
  }),]
  
  # step 3b: estimate principal curve on the meta data
  set.seed(10)
  prin_res <- princurve::principal_curve(x = mat_meta, 
                                         start = mat_meta)
  prin_proj <- princurve::project_to_curve(mat, 
                                           s = prin_res$s)
  
  tmp <- mat[prin_proj$ord,]
  
  data.frame(cell = rownames(tmp), ord = 1:nrow(tmp))
}