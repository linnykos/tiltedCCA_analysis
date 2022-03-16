.compute_fate_direction <- function(mat_firsttime,
                                    mat_secondtime,
                                    pseudotime_rank,
                                    snn_mat,
                                    verbose = 0){
  stopifnot(length(pseudotime_rank) == nrow(snn_mat),
            nrow(snn_mat) == ncol(snn_mat),
            length(pseudotime_rank) == nrow(mat_firsttime),
            length(pseudotime_rank) == nrow(mat_secondtime),
            ncol(mat_firsttime) == ncol(mat_secondtime))
  
  cell_idx <- which(!is.na(pseudotime_rank))
  fit_list <- lapply(1:length(cell_idx), function(k){
    if(verbose > 0 && k %% floor(length(cell_idx)/10) == 0) cat('*')
    
    i <- cell_idx[k]
    nn_idx <- tiltedCCA:::.nonzero_col(mat = snn_mat, col_idx = i, bool_value = F)
    nn_idx <- intersect(nn_idx, cell_idx)
    nn_idx <- unique(c(nn_idx, i))
    cell_rank <- pseudotime_rank[i]
    nn_rank <- pseudotime_rank[nn_idx]
    nn_val <- nn_rank - cell_rank
    nn_val <- nn_val/max(abs(nn_val))
    
    fit_vec <- sapply(nn_idx, function(j){
      df <- data.frame(first = mat_firsttime[i,],
                       second = mat_secondtime[j,])
      lm_res <- stats::lm(second ~ first, data = df)
      summary(lm_res)$r.squared
    })
    names(fit_vec) <- as.character(nn_idx)
      
    fate_rank <- nn_val[which.max(fit_vec)]
    r2_fit <- max(fit_vec)
    r2_original <- fit_vec[as.character(i)]
    list(fate_rank = fate_rank,
         r2_fit = r2_fit,
         r2_original = r2_original)
  })
  
  fate_vec_full <- rep(NA, n)
  fate_vec_full[cell_idx] <- sapply(fit_list, function(x){x$fate_rank})
  r2_vec_full <- rep(NA, n)
  r2_vec_full[cell_idx] <- sapply(fit_list, function(x){x$r2_fit})
  r2_org_full <- rep(NA, n)
  r2_org_full[cell_idx] <- sapply(fit_list, function(x){x$r2_original})
  
  list(fate_vec = fate_vec_full,
       r2_vec = r2_vec_full,
       r2_org = r2_org_full)
}
