form_snn_graph <- function(obj, 
                           clustering, 
                           selected_clusters,
                           k = 5){
  # extract relevant embeddings, both for RNA and ATAC
  tmp <- multiomicCCA:::.prepare_embeddings(obj, 
                                            data_1 = T, 
                                            data_2 = F, 
                                            add_noise = F, 
                                            center = F, 
                                            renormalize = F)
  rna_embedding <- tmp[["everything"]]
  rna_embedding <-.form_metacell_matrix(rna_embedding, clustering)
  l2_vec <- apply(rna_embedding, 1, multiomicCCA:::.l2norm)
  rna_embedding <- multiomicCCA:::.mult_vec_mat(1/l2_vec, rna_embedding)
  rna_embedding <- rna_embedding[rownames(rna_embedding) %in% selected_clusters,]
  
  tmp <- multiomicCCA:::.prepare_embeddings(obj, 
                                            data_1 = F, 
                                            data_2 = T, 
                                            add_noise = F, 
                                            center = F, 
                                            renormalize = F)
  atac_embedding <- tmp[["everything"]]
  atac_embedding <-.form_metacell_matrix(atac_embedding, clustering)
  center_vec <- matrixStats::colMeans2(atac_embedding)
  sd_vec <- matrixStats::colSds(atac_embedding)
  for(j in 1:ncol(atac_embedding)){
    atac_embedding[,j] <- (atac_embedding[,j]-center_vec[j])/sd_vec[j]
  }
  l2_vec <- apply(atac_embedding, 1, multiomicCCA:::.l2norm)
  atac_embedding <- multiomicCCA:::.mult_vec_mat(1/l2_vec, atac_embedding)
  atac_embedding <- atac_embedding[rownames(atac_embedding) %in% selected_clusters,]
  
  # compute cosine distances among the k
  rna_dist <- acos(lsa::cosine(t(rna_embedding)))/pi
  atac_dist <- acos(lsa::cosine(t(rna_embedding)))/pi
  dist_mat <- sqrt(rna_dist^2 + atac_dist^2)
  dist_mat2 <- dist_mat
  rownames(dist_mat2) <- NULL
  colnames(dist_mat2) <- NULL
  g <- igraph::graph_from_adjacency_matrix(dist_mat2,
                                           mode = "undirected",
                                           weighted = T)
  g2 <- igraph::mst(g, algorithm = "prim")
  # ensures connectivity
  mst_edge_mat <- igraph::as_edgelist(g2)
  
  # while loop
  adj_mat <- .form_snn_from_edgelists(mst_edge_mat, 
                                      dist_mat,
                                      k = k,
                                      rowname_vec = rownames(dist_mat))
  
  n <- nrow(adj_mat)
  snn <- adj_mat
  for(i in 1:n){
    idx <- which(adj_mat[i,] != 0)
    snn[i,idx] <- dist_mat[i,idx]
  }
  
  # output both snn and adj_mat
  list(snn = snn, adj_mat = adj_mat, dist_mat = dist_mat)
}

form_transition <- function(snn, 
                            lazy_param = 0.85,
                            teleport_param = 0.99){
  n <- nrow(snn)
  P <- snn
  
  for(i in 1:n){
    idx <- which(P[i,] != 0)
    min_val <- min(P[i,idx])
    max_val <- max(P[i,idx])
    P[i,idx] <- exp(-(P[i,idx]- min_val)/max_val)
  }
  P <- diag(1/matrixStats::rowSums2(P)) %*% P
  
  P <- lazy_param*P + (1-lazy_param)*diag(ncol(P))
  P <- teleport_param*P + (1-teleport_param)*matrix(1/n, n, n)
  
  rownames(P) <- rownames(snn)
  colnames(P) <- colnames(snn)
  
  P
}

# try the diffusion distance, see 
# https://www.stat.berkeley.edu/~mmahoney/s15-stat260-cs294/Lectures/lecture15-12mar15.pdf
# https://academic.oup.com/bioinformatics/article/31/18/2989/241305
# https://arxiv.org/pdf/1406.0013.pdf
# https://arxiv.org/pdf/0811.0121.pdf
# https://mathworld.wolfram.com/LeftEigenvector.html
# https://mathworld.wolfram.com/Eigenvector.html
extract_eigen <- function(P, dims = 1:ncol(P), check = F){
  right_eigen <- eigen(P, symmetric = F)
  left_eigen <- eigen(t(P), symmetric = F)
  if(any(Re(left_eigen$vectors[,1]) < 0)) left_eigen$vectors[,1] <- -left_eigen$vectors[,1]
  eigenvalues <- Re(left_eigen$values)
  
  if(check){
    n <- nrow(P)
    stopifnot(sum(abs(P%*%right_eigen$vectors[,2] - right_eigen$values[2]*right_eigen$vectors[,2])) <= 1e-6)
    stopifnot(sum(abs(left_eigen$vectors[,2]%*%P - left_eigen$values[2]*left_eigen$vectors[,2])) <= 1e-6)
    stopifnot(sum(abs(right_eigen$vectors[,1]+rep(1/sqrt(n)))) <= 1e-6 |
                sum(abs(right_eigen$vectors[,1]-rep(1/sqrt(n)))) <= 1e-6)
    stopifnot(sum(abs(sort(right_eigen$values) - sort(left_eigen$values))) <= 1e-6)
    
    stopifnot(abs(sum(left_eigen$vectors[,2]^2) - 1) <= 1e-6)
  }
  
  # normalize
  left_vector <- Re(left_eigen$vectors)
  for(i in 2:ncol(left_vector)){
    left_vector[,i] <- left_vector[,i]*sqrt(left_vector[,1])
  }
  right_vector <- Re(right_eigen$vectors)
  for(i in 1:ncol(right_vector)){
    right_vector[,i] <- right_vector[,i]/sqrt(left_vector[,1])
  }
  
  list(eigenvalues = eigenvalues[dims], 
       left_vector = left_vector[,dims,drop = F],
       right_vector = right_vector[,dims,drop = F])
}

# [[for debugging only]]
diffusion_distance_singleton <- function(eigenvalues, right_vector, 
                                         idx1, idx2, 
                                         time_vec = 1:min(40, length(eigenvalues))){
  sqrt(sum(sapply(time_vec, function(x){
    sum((c(eigenvalues[-1])^x*(right_vector[idx1,-1] - right_vector[idx2,-1]))^2)
  })))
}

diffusion_distance <- function(eigenvalues, 
                               right_vector, 
                               time_vec = 1:min(40, length(eigenvalues))){
  time_list <- lapply(time_vec, function(x){
    multiomicCCA:::.mult_mat_vec(right_vector[,-1,drop = F], eigenvalues[-1]^x)
  })
  
  identity_list <- lapply(time_list, function(mat){
    n <- nrow(mat)
    tmp <- apply(mat, 1, function(y){multiomicCCA:::.l2norm(y)^2})
    matrix(rep(tmp, each = n), ncol = n, nrow = n)
  })
  
  crossprod_list <- lapply(time_list, function(mat){
    tcrossprod(mat)
  })
  
  dist_list <- lapply(1:length(time_vec), function(i){
    identity_list[[i]] + t(identity_list[[i]]) - 2*crossprod_list[[i]]
  })
  
  dist_mat <- Reduce('+', dist_list)
  diag(dist_mat) <- 0
  dist_mat <- sqrt(dist_mat)
  dist_mat
}

###################################

.form_metacell_matrix <- function(dat, clustering, func = median){
  stopifnot(is.character(clustering), length(clustering) == nrow(dat))
  
  uniq_clust <- sort(unique(clustering))
  clust_mat <- t(sapply(uniq_clust, function(clust){
    idx <- which(clustering == clust)
    apply(dat[idx,,drop = F], 2, func)
  }))
  rownames(clust_mat) <- uniq_clust
  
  clust_mat
}


.form_snn_from_edgelists <- function(mst_edge_mat, 
                                     dist_mat,
                                     k,
                                     rowname_vec){
  stopifnot(k <= ncol(dist_mat))
  
  n <- nrow(dist_mat)
  
  adj_mat1 <- matrix(0, n, n)
  for(i in 1:nrow(mst_edge_mat)){
    idx1 <- mst_edge_mat[i,1]
    idx2 <- mst_edge_mat[i,2]
    
    adj_mat1[idx1, idx2] <- 1
    adj_mat1[idx2, idx1] <- 1
  }
  diag(adj_mat1) <- 0
  
  adj_mat2 <- matrix(0, n, n)
  for(i in 1:nrow(dist_mat)){
    adj_mat2[i, order(dist_mat[i,], decreasing = F)[1:(k+1)]] <- 1
  }
  diag(adj_mat2) <- 0
  adj_mat2 <- adj_mat2*t(adj_mat2)
  
  adj_mat <- adj_mat1 + adj_mat2
  adj_mat[adj_mat > 0] <- 1
  
  rownames(adj_mat) <- rowname_vec
  colnames(adj_mat) <- rowname_vec
  
  adj_mat
}