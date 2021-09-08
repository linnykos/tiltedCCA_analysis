form_snn_graph <- function(obj, cell_idx, 
                           k_max = 50,
                           data_1 = T,
                           data_2 = F){
  # extract relevant embeddings, both for RNA and ATAC
  embedding <- multiomicCCA:::.prepare_embeddings(obj, 
                                                  data_1 = data_1, 
                                                  data_2 = data_2, 
                                                  add_noise = F, 
                                                  center = T, 
                                                  renormalize = F)
  
  # tweak said embedding
  mat <- embedding[["everything"]][cell_idx,]
  l2_vec <- apply(mat, 1, multiomicCCA:::.l2norm)
  mat <- multiomicCCA:::.mult_vec_mat(1/l2_vec, mat)
  
  # compute cosine distances among the k_max
  nn_res <- RANN::nn2(mat, k = k_max)
  # convert to angular distances
  nn_res$nn.dist <- 2*acos(1-(nn_res$nn.dist^2)/2)/pi
  
  # compute MST
  # construct graph
  n <- nrow(mat)
  dist_mat <- matrix(0, n, n)
  for(i in 1:n){
    dist_mat[i, nn_res$nn.idx] <- nn_res$nn.dist
  }
  g <- igraph::graph_from_adjacency_matrix(dist_mat,
                                           mode = "undirected",
                                           weighted = T)
  g2 <- igraph::mst(g, algorithm = "prim")
  # ensures connectivity
  mst_edge_mat <- igraph::as_edgelist(g2)
  
  # compute k_fixed SNN
  
  # while loop
  
  # output both snn and adj_mat
}

compute_diffusion <- function(snn){
  
}