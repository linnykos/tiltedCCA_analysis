rm(list=ls())
load("tests/assets/test_data2.RData")
mat_1 <- test_data$mat_1
mat_2 <- test_data$mat_2
seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1))
seurat_obj[["RNA"]]@var.features <- colnames(mat_1)

n <- nrow(mat_1)
large_clustering_1 <- test_data$clustering_1
large_clustering_2 <- test_data$clustering_2
multiSVD_obj <- create_multiSVD(mat_1 = mat_1, mat_2 = mat_2,
                                dims_1 = 1:2, dims_2 = 1:2,
                                center_1 = F, center_2 = F,
                                normalize_row = T,
                                normalize_singular_value = F,
                                recenter_1 = F, recenter_2 = F,
                                rescale_1 = F, rescale_2 = F,
                                scale_1 = F, scale_2 = F)
multiSVD_obj <- form_metacells(input_obj = multiSVD_obj,
                               large_clustering_1 = large_clustering_1, 
                               large_clustering_2 = large_clustering_2,
                               num_metacells = NULL)
multiSVD_obj <- compute_snns(input_obj = multiSVD_obj,
                             latent_k = 2,
                             num_neigh = 10,
                             bool_cosine = T,
                             bool_intersect = T,
                             min_deg = 10)
multiSVD_obj <- tiltedCCA(input_obj = multiSVD_obj,
                          verbose = F)
multiSVD_obj <- tiltedCCA_decomposition(multiSVD_obj)

res <- postprocess_graph_alignment(
  input_obj = multiSVD_obj,
  bool_use_denoised = F,
  bool_include_intercept = T,
  bool_use_metacells = F,
  input_assay = 1,
  return_everything_mat = T,
  seurat_obj = seurat_obj,
  seurat_assay = "RNA",
  seurat_slot = "counts",
  verbose = 0)

everything_mat <- res$everything_mat
colnames(everything_mat)

idx <- which(colnames(everything_mat) %in% c("g16", "g2"))
selected_mat <- everything_mat[,idx]
other_mat <- everything_mat[,-idx]
cor_vec <- sapply(1:ncol(other_mat), function(j){
  .linear_regression(bool_include_intercept = T,
                     bool_center_x = T,
                     bool_center_y = T,
                     bool_scale_x = T,
                     bool_scale_y = T,
                     return_type = "r_squared", 
                     x_mat = selected_mat,
                     y_vec = other_mat[,j])
})

df <- data.frame(cbind(other_mat[,1], selected_mat))
colnames(df)[1] <- "y"
lm_res <- stats::lm(y ~ ., data = df)
plot(df[,1], lm_res$fitted.values)
