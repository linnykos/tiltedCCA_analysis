latent_k = 20
num_neigh = 30
bool_cosine = T
bool_intersect = T
min_deg = 10
verbose = 3
tol = 1e-4
input_obj = multiSVD_obj

param <- tiltedCCA:::.get_param(input_obj)

metacell_clustering <- tiltedCCA:::.get_metacell(input_obj,
                                                 resolution = "cell", 
                                                 type = "list", 
                                                 what = "metacell_clustering")
n <- nrow(tiltedCCA:::.get_SVD(input_obj)$u)
if(!all(is.null(metacell_clustering))){
  averaging_mat <- tiltedCCA:::.generate_averaging_matrix(metacell_clustering = metacell_clustering,
                                                          n = n)
} else {
  averaging_mat <- NULL
}

input_obj <- tiltedCCA:::.set_defaultAssay(input_obj, assay = 1)
dimred_1 <- tiltedCCA:::.get_postDimred(input_obj, averaging_mat = averaging_mat)
input_obj <- tiltedCCA:::.set_defaultAssay(input_obj, assay = 2)
dimred_2 <- tiltedCCA:::.get_postDimred(input_obj, averaging_mat = averaging_mat)

snn_1 <- tiltedCCA:::.form_snn_mat(mat = dimred_1,
                                   num_neigh = num_neigh,
                                   bool_cosine = bool_cosine, 
                                   bool_intersect = bool_intersect, 
                                   min_deg = min_deg,
                                   tol = tol,
                                   verbose = verbose)
snn_2 <- tiltedCCA:::.form_snn_mat(mat = dimred_2,
                                   num_neigh = num_neigh,
                                   bool_cosine = bool_cosine, 
                                   bool_intersect = bool_intersect, 
                                   min_deg = min_deg,
                                   tol = tol,
                                   verbose = verbose)

clustering_1 <- tiltedCCA:::.get_metacell(input_obj,
                                          resolution = "metacell", 
                                          type = "factor",
                                          what = "large_clustering_1")
clustering_2 <- tiltedCCA:::.get_metacell(input_obj,
                                          resolution = "metacell", 
                                          type = "factor",
                                          what = "large_clustering_2")

#############3

n <- nrow(snn_1)
i <- 3014

nn_1 <- tiltedCCA:::.nonzero_col(snn_1, 
                                 col_idx = i,
                                 bool_value = F)
nn_2 <- tiltedCCA:::.nonzero_col(snn_2, 
                                 col_idx = i,
                                 bool_value = F)

prior_1 <- table(clustering_1[nn_2]); prior_1 <- prior_1/sum(prior_1)
prior_2 <- table(clustering_2[nn_1]); prior_2 <- prior_2/sum(prior_2)

idx_all <- sort(unique(c(nn_1, nn_2)))
idx_df <- data.frame(idx = idx_all, 
                     clustering_1 = clustering_1[idx_all],
                     clustering_2 = clustering_2[idx_all])
na_idx <- unique(c(which(is.na(idx_df$clustering_1)), which(is.na(idx_df$clustering_2))))
if(length(na_idx) > 0){
  idx_df <- idx_df[-na_idx,]
}
idx_all <- idx_df$idx
if(length(idx_all) < num_neigh) return(idx_df$idx)

obs_tab <- table(clustering_1[idx_all], clustering_2[idx_all])
obs_tab <- tiltedCCA:::.remove_all_zeros_rowcol(obs_tab)
desired_tab <- tiltedCCA:::.l2_selection_lp(num_neigh = num_neigh,
                                            obs_tab = obs_tab,
                                            prior_1 = prior_1[names(prior_1) %in% rownames(obs_tab)],
                                            prior_2 = prior_2[names(prior_2) %in% colnames(obs_tab)])

tiltedCCA:::.select_l2_cells(desired_tab = desired_tab,
                             idx_df = idx_df)