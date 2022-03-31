.shape_constrained_fit_search <- function(vec, 
                                          grid_vec,
                                          metric){
  fit_list <- lapply(c("concave", "convex", "decreasing", "increasing"), 
                     function(setting){
                       .shape_constrained_fit(vec = vec,
                                              grid_vec = grid_vec,
                                              metric = metric,
                                              setting = setting)
                     })
  mse_vec <- sapply(fit_list, function(fit_vec){
    .error_metric(metric = metric,
                  vec1 = fit_vec,
                  vec2 = vec)
  })
  
  
}

.shape_constrained_fit <- function(vec,
                                   grid_vec, # default: round(seq(0,length(vec),length.out=7))[-c(1,7)]
                                   metric, # default: 2
                                   setting){
  stopifnot(length(setting) == 1, 
            setting %in% c(
              "concave", "concave-convex",
              "convex", "convex-concave",
              "decreasing", "increasing"))
  n <- length(vec)
  if(setting %in% c("decreasing", "increasing")){
    bool_decreasing <- setting == "decreasing"
    res <- UniIsoRegression::reg_1d(vec,
                                    w_vec = rep(1, n),
                                    metric = metric,
                                    unimodal = F,
                                    decreasing = bool_decreasing)
    
  } else if(setting %in% c("concave", "convex")){
    if(setting == "convex") vec <- -vec
    res <- UniIsoRegression::reg_1d(vec,
                                    w_vec = rep(1, n),
                                    metric = metric,
                                    unimodal = T,
                                    decreasing = F)
    if(setting == "convex") res <- -res
    
  } else if(setting %in% c("concave-convex", "convex-concave")){
    combn_mat <- utils::combn(grid_vec, 2)
    mse_list <- lapply(1:ncol(combn_mat), function(j){
      vec_list <- list(vec[1:combn_mat[1,j]],
                       vec[(combn_mat[1,j]+1):combn_mat[2,j]],
                       vec[(combn_mat[2,j]+1):n])
      if(setting == "concave-convex"){
        decreasing_vec <- c(F, T, F)
      } else {
        decreasing_vec <- c(T, F, T)
      }
      
      res <- unlist(lapply(1:3, function(j){
        UniIsoRegression::reg_1d(vec_list[[j]],
                                 w_vec = rep(1, length(vec_list[[j]])),
                                 metric = metric,
                                 unimodal = F,
                                 decreasing = decreasing_vec[j])
      }))
      error_val <- .error_metric(metric, res, vec)
      
      list(fit = res, error = error_val)
    })
    
    idx <- which.min(sapply(mse_list, function(x){x$error}))
    res <- mse_list[[idx]]$fit
  }
  
  res
}

.error_metric <- function(metric, vec1, vec2){
  stopifnot(metric %in% c(1,2,3))
  if(metric == 1){
    res <- sum(abs(vec1 - vec2))
  } else if(metric == 2){
    res <- sqrt(sum((vec1 - vec2)^2))
  } else if(metric == 3){
    res <- max(abs(vec1 - vec2))
  } 
  
  res
}

.substantial_unimodal <- function(vec,
                                  threshold_perc){ # default: 0.1
  n <- length(vec)
  threshold_val <- threshold_perc * n
  sign_vec <- sign(diff(vec))
  positive_region <- range(which(sign_vec > 0))
  negative_region <- range(which(sign_vec < 0))
  
  diff(positive_region) > threshold_val & diff(negative_region) > threshold_val
}
