plot_scores_heatmap <- function(score_mat, membership_vec = NA, 
                                num_col = 10,
                                bool_rownormalize_before = T,
                                bool_center = T, bool_scale = T,
                                bool_rownormalize_after = F,
                                bool_log = F, scaling_power = 1, luminosity = F,
                                ...){
  
  n <- nrow(score_mat)
  if(bool_rownormalize_before){
    rowsum_vec <- matrixStats::rowSums2(score_mat)
    score_mat <- multiomicCCA:::.mult_vec_mat(1/rowsum_vec, score_mat)
  }
  
  if(bool_center | bool_scale){
    score_mat <- scale(score_mat, center = bool_center, scale = bool_scale)
  }
  
  if(bool_rownormalize_after){
    rowsum_vec <- matrixStats::rowSums2(score_mat)
    score_mat <- multiomicCCA:::.mult_vec_mat(1/rowsum_vec, score_mat)
  }
  
  if(bool_log){
    score_mat <- log(abs(score_mat)+1)*sign(score_mat)
  }
  zlim <- range(score_mat)
  
  # construct colors. green is negative
  max_val <- max(abs(zlim)); min_val <- max(min(abs(zlim)), 1e-3)
  col_vec_neg <- .colorRamp_custom(c(0.584, 0.858, 0.564), c(1,1,1), num_col,
                                   luminosity = luminosity)
  break_vec_neg <- -max_val*seq(1, 0, length.out = num_col+2)^scaling_power
  break_vec_neg <- break_vec_neg[-length(break_vec_neg)]
  
  # red is positive
  col_vec_pos <- .colorRamp_custom(c(1,1,1), c(0.803, 0.156, 0.211), num_col,
                                   luminosity = luminosity)
  break_vec_pos <- max_val*seq(0, 1, length.out = num_col+2)^scaling_power
  break_vec_pos <- break_vec_pos[-1]
  
  # combine the two
  break_vec <- c(break_vec_neg, break_vec_pos)
  col_vec <- c(col_vec_neg, "white", col_vec_pos)
  
  if(!all(is.na(membership_vec))){
    stopifnot(is.factor(membership_vec), length(membership_vec) == nrow(score_mat))
    
    membership_vec <- as.numeric(membership_vec) ## convert to integers
    idx <- order(membership_vec, decreasing = F)
    breakpoints <- 1-which(abs(diff(sort(membership_vec, decreasing = F))) >= 1e-6)/n
  } else {
    idx <- 1:n
  }
  
  line_func <- function(){
    if(!all(is.na(membership_vec))){
      for(i in 1:length(breakpoints)){
        graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2.1, col = "white")
        graphics::lines(c(-10, 10), rep(breakpoints[i], 2), lwd = 2, lty = 2)
      }
    }
  }
  
  graphics::image(.rotate(score_mat[idx,,drop = F]), col = col_vec,
                  breaks = break_vec, ...)
  line_func()
  
  invisible()
}


.colorRamp_custom <- function(vec1, vec2, length, luminosity){
  mat <- matrix(0, nrow = length, ncol = 3)
  for(i in 1:3){
    mat[,i] <- seq(vec1[i], vec2[i], length.out = length)
  }
  
  if(luminosity){
    luminosity_vec <- apply(mat, 1, function(x){
      0.2126*x[1] + 0.7152*x[2] + 0.0722*x[3]
    })
    
    target_luminosity <- mean(c(luminosity_vec[1], luminosity_vec[length]))
    
    mat <- t(sapply(1:nrow(mat), function(x){
      factor <- min(c(target_luminosity/luminosity_vec[x], 1/mat[x,]))
      mat[x,] * factor
    }))
  }
  
  apply(mat, 1, function(x){
    grDevices::rgb(x[1], x[2], x[3])
  })
}


.rotate = function(a) { t(a[nrow(a):1,]) }