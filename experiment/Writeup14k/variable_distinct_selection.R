variable_distinct_selection <- function(dcca_res,
                                        significance_vec, #larger is more significant
                                        max_variables = 10,
                                        cor_threshold = 0.9,
                                        verbose = T){
  stopifnot(length(names(significance_vec)) == length(significance_vec),
            all(sort(significance_vec) == sort(rownames(dcca_res$svd_1$v))) ||
              all(sort(significance_vec) == sort(rownames(dcca_res$svd_2$v))))
  
  if(verbose) print("Extracting relevant matrices")
  reference_mat <- dcca_res$common_score
  l2_vec <- apply(reference_mat, 2, .l2norm)
  reference_mat <- .mult_mat_vec(reference_mat, 1/l2_vec)
  
  decomp_res <- dcca_decomposition(dcca_res)
  if(all(sort(variable_significance) == sort(rownames(dcca_res$svd_1$v)))){
    distinct_mat <- decomp_res$distinct_mat_1
  } else {
    distinct_mat <- decomp_res$distinct_mat_2
  }
  l2_vec <- apply(distinct_mat, 2, .l2norm)
  distinct_mat <- .mult_mat_vec(distinct_mat, 1/l2_vec)
  
  n <- nrow(distinct_mat)
  
  selected_variables <- numeric(0)
  while(length(selected_variables) < max_variables){
    if(verbose) print(paste0("On iteration: ", length(selected_variables)+1))
    
    if(verbose) print("Selecting variable")
    cor_vec <- sapply(1:ncol(distinct_mat), function(i){
      .variable_correlation(x_mat = reference_mat, 
                            y_mat = distinct_mat[,i])
    })
    candidate_var <- colnames(distinct_mat)[which(cor_vec <= cor_threshold)]
    if(verbose) print(paste0(length(candidate_var), " eligble variables"))
    if(length(candidate_var) == 0) break()
    idx <- candidate_var[which.max(significance_vec[candidate_var])]
    selected_variables <- c(selected_variables, idx)
    reference_mat <- cbind(reference_mat, distinct_mat[,idx])
    distinct_mat <- distinct_mat[,which(!colnames(distinct_mat) %in% idx),drop = F]
    if(ncol(distinct_mat) == 0) break()
  }
  
  selected_variables
}

.variable_correlation <- function(x_mat, y_vec){
  df <- data.frame(x_mat)
  df <- cbind(df, y_vec)
  colnames(df)[ncol(df)] <- "y"
  
  lm_res <- stats::lm(y ~ ., data = df)
  summary(lm_res)$r.squared
}