input_obj = multiSVD_obj
bool_use_denoised = T
bool_include_intercept = T
bool_use_metacells = F
bool_use_both_modalities = T
cell_idx = cell_idx
cor_threshold = 0.8
num_variables = 50
sd_quantile = 0.75
seurat_obj = greenleaf
seurat_assay_1 = "SCT"
seurat_assay_2 = "customGAct"
seurat_slot = "data"
verbose = 2

res <- tiltedCCA:::.smooth_variable_selection_helper(
  bool_include_intercept = bool_include_intercept,
  bool_use_denoised = bool_use_denoised,
  bool_use_metacells = bool_use_metacells,
  cell_idx = cell_idx,
  input_assay = 1,
  input_obj = input_obj,
  sd_quantile = sd_quantile,
  seurat_assay = seurat_assay_1,
  seurat_obj = seurat_obj,
  seurat_slot = seurat_slot,
  verbose = verbose
)
alignment_vec_1 <- res$alignment_vec
everything_mat_1 <- res$everything_mat
sd_vec_1 <- res$sd_vec

res <- tiltedCCA:::.smooth_variable_selection_helper(
  bool_include_intercept = bool_include_intercept,
  bool_use_denoised = bool_use_denoised,
  bool_use_metacells = bool_use_metacells,
  cell_idx = cell_idx,
  input_assay = 2,
  input_obj = input_obj,
  sd_quantile = sd_quantile,
  seurat_assay = seurat_assay_2,
  seurat_obj = seurat_obj,
  seurat_slot = seurat_slot,
  verbose = verbose
)
alignment_vec_2 <- res$alignment_vec
everything_mat_2 <- res$everything_mat
sd_vec_2 <- res$sd_vec
