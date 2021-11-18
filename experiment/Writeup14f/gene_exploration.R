compute_variable_summary <- function(mat, common_mat, 
                                     metacell_clustering,
                                     verbose = 1){
  n <- nrow(mat); p <- ncol(mat)
  
  summary_mat <- t(sapply(1:p, function(j){
    if(verbose == 1 & p > 10 & j %% floor(p/10) == 0) cat('*')
    if(verbose >= 2) print(j)
    
    r_val <- .gene_regression(x_vec = common_mat[,j],
                              y_vec = mat[,j])
    kl_val <- .separation_value(vec = mat[,j],
                                metacell_clustering = metacell_clustering)
    
    c(r_val = r_val, kl_val = kl_val)
  }))
 
  colnames(summary_mat) <- c("r_squared", "kl_div")
  rownames(summary_mat) <- colnames(mat)
  
  summary_mat
}

.gene_regression <- function(x_vec, y_vec){
  x_vec <- scale(x_vec)
  y_vec <- scale(y_vec)
  
  df <- data.frame(x = x_vec, y = y_vec)
  
  lm_res <- stats::lm(y ~ ., data = df)
  summary(lm_res)$r.squared
}

.separation_value <- function(vec, metacell_clustering){
  clust_res <- mclust::Mclust(vec, modelNames = "V", verbose = F)$classification
  clust_res <- factor(clust_res)
  
  baseline_proportion <- table(metacell_clustering)
  baseline_proportion <- baseline_proportion/sum(baseline_proportion)
  
  kl_vec <- sapply(levels(clust_res), function(clust_val){
    idx <- which(clust_res == clust_val)
    clust_proportion <- table(metacell_clustering[idx])
    clust_proportion <- clust_proportion/sum(clust_proportion)
    
    multiomicCCA:::.kl_divergence(query_dist = clust_proportion,
                                  reference_dist = baseline_proportion)
  })
  
  mean(kl_vec)
}