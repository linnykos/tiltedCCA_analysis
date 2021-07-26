k_max <- 50
# c_eig <- multiomicCCA::compute_laplacian(rna_frnn$c_g, k_max = k_max, rowname_vec = row_vec, 
#                                         colname_vec = paste0("clap_", 1:k_max), verbose = F)

###############

g_mat <- rna_frnn$c_g
n <- nrow(g_mat)
g_mat <- multiomicCCA:::.symmetrize_sparse(g_mat, set_ones = T)

deg_vec <- sparseMatrixStats::rowSums2(g_mat)
invdeg_mat <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1/deg_vec)
g_mat <- invdeg_mat %*% g_mat %*% invdeg_mat