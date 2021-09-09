detect_changepoints <- function(obj, 
                                prin_df,
                                threshold = 0.01){
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
  solution_list <- lapply(1:p, function(j){
    if(p > 10 && j %% floor(p/10) == 0) cat('*')
    # find changepoints
    res <- breakfast::sol.idetect(mat[odd_idx,j])
    if(length(res$solution.path) == 0) return(NA)
    
    changepoints <- sort(res$solution.path[1:min(2,length(res$solution.path))])
    
    # test changepoints
    .wilcox_tests(mat[even_idx,j], 
                  changepoints = changepoints,
                  threshold = threshold)
  })
  names(solution_list) <- colnames(mat)
  
  solution_list <- solution_list[sapply(solution_list, function(x){!all(is.na(x))})]
  solution_list
}

####################3

.wilcox_tests <- function(vec, 
                          changepoints, 
                          threshold = 0.01){
  stopifnot(length(changepoints) <= 2)
  n <- length(vec)
  
  if(length(changepoints) == 1){
    res1 <- stats::wilcox.test(x = vec[1:changepoints],
                               y = vec[changepoints:n])
    if(res1$p.value <= threshold){
      mean1 <- mean(vec[1:changepoints])
      mean2 <- mean(vec[changepoints:n])
      
      if(mean1 > mean2){
        return(list(start = 1,
                    end = changepoint,
                    diff = mean1 - mean2,
                    pval = res1$p.value))
      } else {
        return(list(start = changepoint,
                    end = n,
                    diff = mean2 - mean1,
                    pval = res1$p.value))
      }
    } else {
      return(NA)
    }
    
  } else {
    stopifnot(changepoints[1] < changepoints[2])
    res_list <- vector("list", 3)
    # test A vs BC
    res_list[[1]] <- stats::wilcox.test(x = vec[1:changepoints[1]],
                                        y = vec[changepoints[1]:n],
                                        alternative = "two.sided")
    
    # test AB vs C
    res_list[[2]] <- stats::wilcox.test(x = vec[1:changepoints[2]],
                                        y = vec[changepoints[2]:n],
                                        alternative = "two.sided")
    
    # test AC vs B
    res_list[[3]] <- stats::wilcox.test(x = vec[c(1:changepoints[1], changepoints[2]:n)],
                                        y = vec[changepoints[1]:changepoints[2]],
                                        alternative = "less")
    
    pval_vec <- sapply(res_list, function(x){x$p.value})
    if(all(pval_vec > threshold)) return(NA)
    idx <- which.min(pval_vec)
    
    if(idx == 1){
      mean1 <- mean(vec[1:changepoints[1]])
      mean2 <- mean(vec[changepoints[1]:n])
      
      if(mean1 > mean2){
        return(list(start = 1,
                    end = changepoints[1],
                    diff = mean1 - mean2,
                    pval = res_list[[1]]$p.value))
      } else {
        return(list(start = changepoints[1],
                    end = n,
                    diff = mean2 - mean1,
                    pval = res_list[[1]]$p.value))
      }
    } else if(idx == 2){
      mean1 <- mean(vec[1:changepoints[2]])
      mean2 <- mean(vec[changepoints[2]:n])
      
      if(mean1 > mean2){
        return(list(start = 1,
                    end = changepoints[2],
                    diff = mean1 - mean2,
                    pval = res_list[[2]]$p.value))
      } else {
        return(list(start = changepoints[2],
                    end = n,
                    diff = mean2 - mean1,
                    pval = res_list[[2]]$p.value))
      }
    } else{
      mean1 <- mean(vec[c(1:changepoints[1], changepoints[2]:n)])
      mean2 <- mean(vec[changepoints[1]:changepoints[2]])
      
      if(mean1 < mean2){
        return(list(start = changepoints[1],
                    end = changepoints[2],
                    diff = mean2 - mean1,
                    pval = res_list[[3]]$p.value))
      } else {
        return(NA)
      }
    }
  }
}