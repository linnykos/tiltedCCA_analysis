rm(list=ls())
library(Seurat)
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca2.RData")
load("../../../../out/Writeup14i/Writeup14i_citeseq_pbmc224_dcca_variable_full_de.RData")
dcca_res <- dcca_res2
source("../Writeup14f/gene_exploration.R")

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)

summary_mat <- compute_variable_summary(mat = mat_2b, 
                                        common_mat = dcca_decomp$common_mat_2,
                                        metacell_clustering = NA)
colnames(summary_mat) <- c("r_squared", "p_value")
for(i in 1:nrow(summary_mat)){
  gene_name <- rownames(summary_mat)[i]
  
  max_celltype_idx <- max(combn_mat)
  celltype_pval <- sapply(1:max_celltype_idx, function(celltype){
    idx <- which(combn_mat == celltype, arr.ind = T)[,2]
    stats::quantile(sapply(idx, function(j){
      idx <- which(rownames(de_list[[j]]) == gene_name)
      if(length(idx) == 0) return(1)
      de_list[[j]][idx, "p_val"]
    }), probs = 0.75)
  })
  
  summary_mat[i,"p_value"] <- min(celltype_pval)
}
summary_mat[,"p_value"] <- -log10(summary_mat[,"p_value"])

####################################3

.variable_correlation <- function(x_mat, y_vec){
  df <- data.frame(x_mat)
  df <- cbind(df, y_vec)
  colnames(df)[ncol(df)] <- "y"
  
  lm_res <- stats::lm(y ~ ., data = df)
  summary(lm_res)$r.squared
}

variable_significance <- summary_mat[,"p_value"] 
max_variables = 10
num_null = 500
significant_threshold = 0.05
verbose = T

if(verbose) print("Extracting relevant matrices")
reference_mat <- dcca_res$common_score

if(length(variable_significance) == nrow(dcca_res$svd_1$v) &&
   all(sort(variable_significance) == sort(rownames(dcca_res$svd_1$v)))){
  distinct_mat <- tcrossprod(multiomicCCA:::.mult_mat_vec(dcca_res$svd_1$u, dcca_res$svd_1$d), dcca_res$svd_1$v)
} else {
  distinct_mat <- tcrossprod(multiomicCCA:::.mult_mat_vec(dcca_res$svd_2$u, dcca_res$svd_2$d), dcca_res$svd_2$v)
}
colnames(distinct_mat) <- rownames(dcca_res$svd_2$v)
n <- nrow(distinct_mat)

if(verbose) print("Generating the null vectors")
sd_val <- min(apply(reference_mat, 1, stats::sd))
null_mat <- sapply(1:num_null, function(i){
  col_weights <- abs(stats::rnorm(ncol(reference_mat)))
  col_weights <- col_weights/sum(col_weights)
  row_weights <- stats::rnorm(nrow(reference_mat), mean = 1, sd = 1/3)
  multiomicCCA:::.mult_vec_mat(row_weights, reference_mat %*% col_weights)
})

if(verbose) print("Computing threshold")
threshold_vec <- sapply(1:num_null, function(i){
  if(i %% floor(num_null/10) == 0) cat('*')
  .variable_correlation(x_mat = reference_mat, 
                        y_vec = null_mat[,i])
})
print(quantile(threshold_vec))
threshold_val <- stats::quantile(threshold_vec, probs = significant_threshold)
if(verbose) print(paste0("Threshold is R-squared of ", round(threshold_val, 2)))

selected_variables <- numeric(0)
while(length(selected_variables) < max_variables){
  set.seed(length(selected_variables))
  if(verbose) print(paste0("On iteration: ", length(selected_variables)+1))
  
  if(verbose) print("Selecting variable")
  test_vec <- sapply(1:ncol(distinct_mat), function(i){
    if(i %% floor(ncol(distinct_mat)/10) == 0) cat('*')
    .variable_correlation(x_mat = reference_mat, 
                          y_vec = distinct_mat[,i])
  })
  print(quantile(test_vec))
  candidate_var <- colnames(distinct_mat)[which(test_vec <= threshold_val)]
  if(verbose) print(paste0(length(candidate_var), " eligble variables"))
  if(length(candidate_var) == 0) break()
  idx <- candidate_var[which.max(variable_significance[candidate_var])]
  selected_variables <- c(selected_variables, idx)
  reference_mat <- cbind(reference_mat, distinct_mat[,idx])
  distinct_mat <- distinct_mat[,which(!colnames(distinct_mat) %in% idx),drop = F]
  if(ncol(distinct_mat) == 0) break()
}
selected_variables

summary_mat[selected_variables,]
