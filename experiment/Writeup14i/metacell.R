form_subclusters <- function(mat, clustering, target_k, 
                             verbose = T){
  df <- .compute_subcluster_splits(clustering = clustering, 
                                   n = nrow(mat), 
                                   target_k = target_k)
  
  tmp <- lapply(1:length(clustering), function(k){
    if(verbose) print(paste0("On cluster ", k))
    .compute_metacells(k = df$num[k], 
                       mat = mat[clustering[[k]],,drop = F],
                       row_indices = clustering[[k]])
  })
  
  do.call(c, tmp)
}

intersect_metacells <- function(metacell_clustering_1, 
                                metacell_clustering_2,
                                min_size = 5,
                                n){
  stopifnot(n >= max(unlist(metacell_clustering_1)),
            n >= max(unlist(metacell_clustering_2)))
  total_1 <- length(metacell_clustering_1)
  total_2 <- length(metacell_clustering_2)
  
  idx_df <- data.frame(modality_1 = rep(1:total_1, each = total_2),
                       modality_2 = rep(1:total_2, times = total_1))
  idx_df$total_idx <- (idx_df$modality_1 - 1)*total_2 + idx_df$modality_2
  
  idx <- matrix(NA, nrow = n, ncol = 2)
  rownames(idx) <- 1:n
  for(i in 1:total_1){
    idx[metacell_clustering_1[[i]],1] <- i
  }
  for(i in 1:total_2){
    idx[metacell_clustering_2[[i]],2] <- i
  }
  na_idx <- which(apply(idx, 1, function(vec){any(is.na(vec))}))
  if(length(na_idx) > 0) idx <- idx[-na_idx,]
  
  total_idx <- (idx[,1]-1)*total_2 + idx[,2]
  uniq_vec <- sort(unique(total_idx))
  idx_df <- idx_df[idx_df$total_idx %in% uniq_vec,]
  stopifnot(length(idx_df$total_idx) == length(uniq_vec))
  uniq_vec <- idx_df$total_idx

  metacell_clustering <- lapply(uniq_vec, function(i){
    which(total_idx == i)
  })
  names(metacell_clustering) <- as.character(uniq_vec)
  metacell_clustering <- metacell_clustering[sapply(metacell_clustering, length) >= min_size]
  uniq_vec <- as.numeric(names(metacell_clustering))
  
  clustering_hierarchy_1 <- sapply(1:total_1, function(i){
    idx <- which(idx_df$modality_1 == i)
    if(length(idx) > 0){
      which(uniq_vec %in% idx_df$total_idx[idx])
    } else NA
  })
  clustering_hierarchy_2 <- sapply(1:total_2, function(i){
    idx <- which(idx_df$modality_2 == i)
    if(length(idx) > 0){
      which(uniq_vec %in% idx_df$total_idx[idx])
    } else NA
  })
  
  
  list(metacell_clustering = metacell_clustering,
       clustering_hierarchy_1 = clustering_hierarchy_1,
       clustering_hierarchy_2 = clustering_hierarchy_2)
}

#######################

.compute_subcluster_splits <- function(clustering, n, target_k){
  stopifnot(is.list(clustering), 
            sum(sapply(clustering, length)) <= n,
            length(clustering) <= target_k,
            target_k > 0, target_k %% 1 == 0)
  tmp <- unlist(clustering)
  stopifnot(all(tmp %% 1 == 0), table(tmp) == 1, all(tmp > 0))
  
  df <- data.frame(total_size = sapply(clustering, length),
                   size = sapply(clustering, length), 
                   num = rep(1, length(clustering)))
  if(all(is.null(names(clustering)))){
    rownames(df) <- 1:length(clustering)
  } else {
    rownames(df) <- names(clustering)
  }
  
  while(sum(df$num) < target_k){
    idx <- which.max(df$size)
    df$num[idx] <- df$num[idx]+1
    df$size[idx] <- df$total_size[idx]/df$num[idx]
  }
  
  df
}

.compute_metacells <- function(k, mat, row_indices){
  stopifnot(length(row_indices) == nrow(mat))
  
  svd_res <- multiomicCCA:::.svd_truncated(mat, K = min(k, ncol(mat)), 
                                           symmetric = F, 
                                           rescale = F, 
                                           mean_vec = T, 
                                           sd_vec = F, 
                                           K_full_rank = F)
  dimred <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
  kmeans_res <- stats::kmeans(dimred, centers = k)
  metacell_clustering <- lapply(1:k, function(kk){
    row_indices[which(kmeans_res$cluster == kk)]
  })
  
  metacell_clustering
}