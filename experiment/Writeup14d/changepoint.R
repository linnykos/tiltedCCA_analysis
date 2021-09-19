detect_changepoints <- function(obj, 
                                prin_df,
                                scale_factor = 1,
                                min_len = floor(nrow(prin_df))/2,
                                thres_seq = c(1.15, 1.5, 2, 3)){
  if(is.list(obj)){
    # step 1: reconstruct the smoothed gene expression values
    mat <- tcrossprod(multiomicCCA:::.mult_mat_vec(obj$svd_1$u, obj$svd_1$d), obj$svd_1$v)
    rownames(mat) <- rownames(obj$svd_1$u)
    colnames(mat) <- rownames(obj$svd_1$v)
  } else {
    stopifnot(is.matrix(obj))
    mat <- obj
  }
  prin_df <- prin_df[order(prin_df$ord, decreasing = F),]
  cell_ord <- sapply(prin_df$cell, function(cell_name){
    which(rownames(mat) == cell_name)
  })
  mat <- mat[cell_ord,,drop = F]
  odd_idx <- seq(1, nrow(mat), by = 2)
  even_idx <- seq(2, nrow(mat), by = 2)
  
  p <- ncol(mat)
  n <- floor(nrow(mat)/2)
  solution_list <- lapply(1:p, function(j){
    if(p > 10 && j %% floor(p/10) == 0) cat('*')
    # find changepoints
    
    k <- 1
    while(TRUE){
      cp_fit <- breakfast::sol.idetect_seq(mat[odd_idx,j])
      model_fit <- breakfast::model.thresh(cp_fit, th_const = thres_seq[k])
      cp_vec <- model_fit$cpts
      if(length(cp_vec) == 0 || all(cp_vec == 0)) return(NA)
      if(any(diff(c(0, cp_vec, n)) > min_len)) break()
      
      k <- k+1
      if(k > length(thres_seq)) return(NA)
    }
    
    omit_idx <- which(diff(c(0, cp_vec, n)) > min_len)
    if(length(omit_idx) > 1){
      tmp <- c(0, cp_vec, n)
      idx <- which(diff(tmp) > min_len)
      omit_idx <- idx[which.min(sapply(idx, function(x){
        mean(mat[even_idx,j][max(tmp[idx],1):tmp[idx+1]])
      }))]
    }
    res <- .construct_intervals(mat[even_idx,j], 
                                cp_vec = cp_vec, 
                                omit_idx = omit_idx)
    
    .assess_difference(mat[even_idx,j], 
                       interval_mat = res$interval_mat,
                       baseline_mean = res$baseline_mean,
                       baseline_sd = res$baseline_sd,
                       scale_factor = scale_factor)
  })
  names(solution_list) <- colnames(mat)
  
  solution_list <- solution_list[sapply(solution_list, function(x){!all(is.na(x))})]
  solution_list
}

####################

.construct_intervals <- function(vec, cp_vec, omit_idx){
  n <- length(vec)
  cp_vec2 <- c(0, cp_vec, n)
  interval_mat <- cbind(cp_vec2[-length(cp_vec2)], cp_vec2[-1])
  omit_interval <- interval_mat[omit_idx,]
  interval_mat <- interval_mat[-omit_idx,,drop = F]
  baseline_mean <- mean(vec[max(omit_interval[1],1):omit_interval[2]])
  if(diff(omit_interval) > 3){
    baseline_sd <- stats:: sd(vec[max(omit_interval[1],1):omit_interval[2]])
  } else {
    baseline_sd <- 0
  }
  colnames(interval_mat) <- c("start", "end")
  
  list(interval_mat = interval_mat, 
       baseline_mean = baseline_mean, 
       baseline_sd = baseline_sd)
}

.assess_difference <- function(vec, 
                               interval_mat,
                               baseline_mean,
                               baseline_sd,
                               scale_factor){
  mean_vec <- sapply(1:nrow(interval_mat), function(i){
    mean(vec[max(interval_mat[i,1],1):interval_mat[i,2]])
  })
  idx <- which(mean_vec > baseline_mean+scale_factor*baseline_sd)
  if(length(idx) > 0) {
    interval_mat <- interval_mat[idx,,drop = F]
    diff_val <- sapply(1:nrow(interval_mat), function(i){
      (mean(vec[max(interval_mat[i,1],1):interval_mat[i,2]]) - baseline_mean)/baseline_sd
    })
    interval_mat <- cbind(interval_mat, diff_val)
    colnames(interval_mat) <- c("start", "end", "diff")
    
    return(interval_mat)
  } else {
    return(NA)
  }
}