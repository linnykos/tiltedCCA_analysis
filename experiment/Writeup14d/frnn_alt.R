construct_frnn_alt <- function(obj, 
                               membership_vec, 
                               nn = 30, 
                               data_1 = T, 
                               data_2 = F,
                               max_subsample_frnn = nrow(obj$common_score),
                               frnn_approx = 0, 
                               radius_quantile = 0.5,
                               bool_matrix = T, 
                               center = T,
                               renormalize = F, 
                               symmetrize = F,
                               verbose = T,
                               normalization_type = "none"){
  stopifnot(frnn_approx >= 0, frnn_approx <= 1,
            length(membership_vec) == nrow(obj$common_score),
            is.factor(membership_vec))
  
  embedding <- multiomicCCA:::.prepare_embeddings(obj, 
                                                  data_1 = data_1, 
                                                  data_2 = data_2, 
                                                  add_noise = F, 
                                                  center = center, 
                                                  renormalize = F)
  
  if(normalization_type == "everything"){
    l2_vec <- apply(embedding[["everything"]], 1, multiomicCCA:::.l2norm)
    for(i in 1:3){
      embedding[[i]] <- multiomicCCA:::.mult_vec_mat(1/l2_vec, embedding[[i]])
    }
  } else if(normalization_type == "itself"){
    for(i in 1:3){
      l2_vec <- apply(embedding[[i]], 1, multiomicCCA:::.l2norm)
      embedding[[i]] <- multiomicCCA:::.mult_vec_mat(1/l2_vec, embedding[[i]])
    }
  } else if(normalization_type == "none") {
    # do nothing
  }
  
  n <- nrow(embedding[[1]])
  
  # construct subsamples
  cell_subidx <- multiomicCCA:::.construct_celltype_subsample(membership_vec, max_subsample_frnn)
  if(length(cell_subidx) < n) {
    membership_vec <- membership_vec[cell_subidx]
  }
  for(i in 1:3){
    embedding[[i]] <- embedding[[i]][cell_subidx,,drop = F]
  }
  n <- nrow(embedding[[1]])
  
  # compute the radius
  vec_print <- c("common", "distinct", "everything")
  vec_rad <- sapply(1:3, function(i){
    if(verbose) print(paste0(Sys.time(),": cLISI: Computing radius -- ", vec_print[i]))
    multiomicCCA:::.compute_radius(embedding[[i]], nn, radius_quantile)
  })
  vec_rad_org <- vec_rad
  names(vec_rad_org) <- c("common", "distinct", "everything")
  vec_rad[1:2] <- max(vec_rad[1:2])
  
  # construct frnn
  list_g <- lapply(1:3, function(i){
    if(verbose) print(paste0(Sys.time(),": cLISI: Construct graph -- ", vec_print[i]))
    multiomicCCA:::.construct_frnn(embedding[[i]], 
                                   radius = vec_rad[i], 
                                   nn = nn, 
                                   frnn_approx = frnn_approx, 
                                   verbose = verbose)
  }) 
  
  for(i in 1:3){
    list_g[[i]] <- multiomicCCA:::.nnlist_to_matrix(list_g[[i]])
    if(symmetrize){
      list_g[[i]] <- multiomicCCA:::.symmetrize_sparse(list_g[[i]], set_ones = F)
    }
    
    # convert back to list form if needed
    if(bool_matrix){
      if(length(rownames(obj$common_score)) != 0){
        rownames(list_g[[i]]) <- rownames(obj$common_score)
        colnames(list_g[[i]]) <- rownames(obj$common_score)
      }
    } else {
      # [[note to self: add a test to make sure this conversion is bijective]]
      list_g[[i]] <- multiomicCCA:::.matrix_to_nnlist(list_g[[i]])
    }
  }
  
  # [[note to self: I'm not sure if we need to output e_g]]
  return(list(c_g = list_g[[1]], d_g = list_g[[2]], e_g = list_g[[3]], 
              membership_vec = membership_vec,
              original_radius = vec_rad_org))
}

#######################

combine_frnn_alt <- function(dcca_obj, 
                             g_1, g_2, 
                             nn, 
                             common_1 = T, 
                             common_2 = T, 
                             keep_nn = T,
                             sampling_type = "adaptive_gaussian",
                             center = T, 
                             renormalize = F,
                             verbose = 0){
  stopifnot(all(dim(g_1) == dim(g_2)))
  
  # extract the relevant embeddings from dcca_obj
  if(verbose > 0) print(paste0(Sys.time(),": Preparing first embedding"))
  embedding_1 <- multiomicCCA:::.prepare_embeddings(dcca_obj, data_1 = T, data_2 = F, 
                                                    add_noise = F, center = center, 
                                                    renormalize = renormalize)
  if(common_1){
    embedding_1 <- embedding_1$common
  } else {
    embedding_1 <- embedding_1$distinct
  }
  
  if(verbose > 0) print(paste0(Sys.time(),": Preparing second embedding"))
  embedding_2 <- multiomicCCA:::.prepare_embeddings(dcca_obj, data_1 = F, data_2 = T, 
                                     add_noise = F, center = center, 
                                     renormalize = renormalize)
  if(common_2){
    embedding_2 <- embedding_2$common 
  } else {
    embedding_2 <- embedding_2$distinct
  }
  
  # symmetrize g_1 and g_2
  if(verbose > 0) print(paste0(Sys.time(),": Symmetrizing"))
  g_1 <- multiomicCCA:::.symmetrize_sparse(g_1, set_ones = F)
  g_2 <- multiomicCCA:::.symmetrize_sparse(g_2, set_ones = F)
  
  # prepare
  if(verbose > 0) print(paste0(Sys.time(),": Converting matrices to list"))
  n <- nrow(g_1)
  nn_idx_1 <- lapply(1:n, function(j){multiomicCCA:::.nonzero_col(g_1, j, bool_value = F)})
  nn_dist_1 <- lapply(1:n, function(j){multiomicCCA:::.nonzero_col(g_1, j, bool_value = T)})
  nn_idx_2 <- lapply(1:n, function(j){multiomicCCA:::.nonzero_col(g_2, j, bool_value = F)})
  nn_dist_2 <- lapply(1:n, function(j){multiomicCCA:::.nonzero_col(g_2, j, bool_value = T)})
  
  # apply the following procedure for each cell n
  nn_idx_all <- vector("list", n); nn_dist_all <- vector("list", n)
  for(i in 1:n){
    if(verbose == 2) {
      print(i)
    } else if(verbose == 1 && n > 10 && i %% floor(n/10) == 0) {
      cat('*')
    }
    
    # intersect
    idx_all <- unique(c(nn_idx_1[[i]], nn_idx_2[[i]]))
    idx_intersect <- intersect(nn_idx_1[[i]], nn_idx_2[[i]])
    if(length(idx_intersect) < nn){
      idx_subset <- setdiff(idx_all, idx_intersect)
      idx_intersect <- c(idx_intersect, sample(idx_subset, size = nn - length(idx_intersect)))
    }
    
    # start tabulating nn_dist
    nn_idx_all[[i]] <- idx_intersect
    nn_dist_all[[i]] <- .compute_distance_from_idx_alt(embedding_1, embedding_2, 
                                                   nn_idx_1[[i]], nn_idx_2[[i]],
                                                   nn_dist_1[[i]], nn_dist_2[[i]],
                                                   start_idx = i, end_idx_vec = idx_intersect)
    # sample these entries
    tmp <- multiomicCCA:::.embedding_resampling(nn_idx_all, nn_dist_all, nn = nn, 
                                 sampling_type = sampling_type, keep_nn = F)
    nn_idx_all <- tmp$id; nn_dist_all <- tmp$dist
    
    # find the nn's
    if(keep_nn){
      if(length(nn_idx_1[[i]]) < nn){
        order_1 <- nn_idx_1[[i]]
      } else {
        order_1 <- order(nn_dist_1[[i]], decreasing = F)[1:nn]
      }
      if(length(nn_idx_2[[i]]) < nn){
        order_2 <- nn_idx_2[[i]]
      } else {
        order_2 <- order(nn_dist_2[[i]], decreasing = F)[1:nn]
      }
      
      tmp_idx <- unique(c(nn_idx_1[[i]][order_1], nn_idx_2[[i]][order_2]))
      tmp_dist <- .compute_distance_from_idx_alt(embedding_1, embedding_2, 
                                             nn_idx_1[[i]], nn_idx_2[[i]],
                                             nn_dist_1[[i]], nn_dist_2[[i]],
                                             start_idx = i, end_idx_vec = tmp_idx)
      zz <- intersect(order(tmp_dist, decreasing = F)[1:nn], which(!tmp_idx %in% nn_idx_all[[i]]))
      nn_idx_all[[i]] <- c(nn_idx_all[[i]], tmp_idx[zz])
      nn_dist_all[[i]] <- c(nn_dist_all[[i]], tmp_dist[zz])
    }
  }
  
  tmp_list <- list(id = nn_idx_all, dist = nn_dist_all)
  res <- multiomicCCA:::.nnlist_to_matrix(tmp_list)
  res <- multiomicCCA:::.symmetrize_sparse(res, set_ones = F)
  
  if(length(rownames(dcca_obj$common_score)) != 0){
    rownames(res) <- rownames(dcca_obj$common_score)
    colnames(res) <- rownames(dcca_obj$common_score)
  }
  
  res
}

# based on the fact that sqrt(d_1^2(x,y) + d_2^2(x,y)) is a distance metric itself
.compute_distance_from_idx_alt <- function(embedding_1, embedding_2, 
                                       nn_idx_1_vec, nn_idx_2_vec,
                                       nn_dist_1_vec, nn_dist_2_vec,
                                       start_idx, end_idx_vec){
  stopifnot(!start_idx %in% end_idx_vec)
  
  res <- sapply(end_idx_vec, function(j){
    tmp <- 0
    z1 <- which(nn_idx_1_vec == j)
    if(length(z1) == 1) {
      tmp <- tmp + nn_dist_1_vec[z1]^2
    } else {
      tmp <- tmp + multiomicCCA:::.l2norm(embedding_1[start_idx,] - embedding_1[j,])^2
    }
    
    z2 <- which(nn_idx_2_vec == j)
    if(length(z2) == 1) {
      tmp <- tmp + nn_dist_2_vec[z2]^2
    } else {
      tmp <- tmp + multiomicCCA:::.l2norm(embedding_2[start_idx,] - embedding_2[j,])^2
    }
    
    sqrt(tmp)
  })
  
  res
}