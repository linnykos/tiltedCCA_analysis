multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(multiSVD_obj, 
                                                    bool_modality_1_full = F,
                                                    bool_modality_2_full = F,
                                                    verbose = F)
reference_dimred <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_dimred")

everything_mat <- t(as.matrix(pbmc[["ADT"]]@scale.data))

min_subsample_cell <- 5000
seurat_obj <- pbmc
seurat_celltype_variable <- "celltype.l2"
verbose = 1
if(!is.null(min_subsample_cell)){
  if(verbose > 1) print("Reducing the number of cells")
  stopifnot(seurat_celltype_variable %in% colnames(seurat_obj@meta.data))
  membership_vec <- seurat_obj@meta.data[,seurat_celltype_variable]
  idx <- tiltedCCA:::construct_celltype_subsample(membership_vec, min_subsample_cell = min_subsample_cell)
  
  ncell_before <- nrow(reference_dimred)
  reference_dimred <- reference_dimred[idx,]
  everything_mat <- everything_mat[idx,]
  ncell_after <- nrow(reference_dimred)
  if(verbose > 1) print(paste0("Reduced the number of cells from ", ncell_before, " to ", ncell_after))
}

cor_vec <- sapply(1:ncol(everything_mat), function(j){
  if(verbose > 1 && ncol(everything_mat) > 10 & j %% floor(ncol(everything_mat)/10) == 0) cat('*')
  if(verbose > 2) print(paste0("Working on variable ", j, " of ", ncol(everything_mat)))
  tiltedCCA:::.linear_regression(bool_include_intercept = T,
                                 bool_center_x = T,
                                 bool_center_y = T,
                                 bool_scale_x = T,
                                 bool_scale_y = T,
                                 return_type = "r_squared", 
                                 x_mat = reference_dimred,
                                 y_vec = everything_mat[,j])
})
names(cor_vec) <- colnames(everything_mat)

round(cor_vec[order(names(cor_vec))],2)
