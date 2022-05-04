rm(list=ls())
load("../../../out/main/citeseq_pbmc224_differential.RData")
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")
source("pbmc_224antibody_colorPalette.R")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# gather all the relevant ingredients
Seurat::DefaultAssay(pbmc) <- "ADT"
logpval_vec <- sapply(rownames(pbmc), function(gene){
  val <- stats::quantile(
    sapply(1:length(adt_de_list), function(j){
      idx <- which(rownames(adt_de_list[[j]]) == gene)
      if(length(idx) == 0) return(1)
      adt_de_list[[j]][idx, "p_val"]
    }), probs = 0.75
  )
  
  -log10(val)
})
# tmp <- sapply(1:length(adt_de_list), function(j){
#   adt_de_list[[j]][which(rownames(adt_de_list[[j]]) == "CD52"), "p_val"]
# })
names(logpval_vec) <- rownames(pbmc)
quantile(logpval_vec)
max_val <- max(logpval_vec[which(!is.infinite(logpval_vec))])
logpval_vec <- pmin(logpval_vec, 300)

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(multiSVD_obj, 
                                                    bool_modality_1_full = F,
                                                    bool_modality_2_full = F,
                                                    verbose = 0)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
reference_dimred <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_dimred")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
common_mat <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_mat")
distinct_mat <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "distinct_mat")
adt_mat <- common_mat + distinct_mat

rsquare_vec <- sapply(1:ncol(adt_mat), function(j){
  tiltedCCA:::.linear_regression(bool_include_intercept = T,
                                 bool_center_x = T,
                                 bool_center_y = T,
                                 bool_scale_x = T,
                                 bool_scale_y = T,
                                 return_type = "r_squared", 
                                 x_mat = reference_dimred,
                                 y_vec = adt_mat[,j])
})
names(rsquare_vec) <- colnames(adt_mat)

######################

# First panel: Select the 10 antibodies with the highest p-vales
panel_1 <- names(logpval_vec)[order(logpval_vec, decreasing = T)[1:10]]

# Second panel: Select the 10 antibodies deviating the most from the common space
panel_2 <- names(rsquare_vec)[order(rsquare_vec, decreasing = F)[1:10]]

# Third panel: Select the 10 antibodies using the proposed strategy, but only consider the deviation from the RNA subspace
verbose <- 1
reference_dimred <- tiltedCCA:::.mult_mat_vec(multiSVD_obj$svd_1$u, multiSVD_obj$svd_1$d)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
common_mat <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_mat")
distinct_mat <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "distinct_mat")
adt_mat <- common_mat + distinct_mat
logpval_vec <- logpval_vec[colnames(adt_mat)]

seurat_obj <- pbmc
seurat_celltype_variable <- "celltype.l2"
membership_vec <- seurat_obj@meta.data[,seurat_celltype_variable]
membership_vec <- membership_vec[which(rownames(seurat_obj@meta.data) %in% rownames(reference_dimred))]
idx <- tiltedCCA:::construct_celltype_subsample(membership_vec, 
                                                min_subsample_cell = 5000)
reference_dimred <- reference_dimred[idx,]
adt_mat <- adt_mat[idx,]

cor_threshold <- 0.85
max_variables <- 10
n <- nrow(adt_mat)
candidate_list <- vector("list", length = max_variables)
selected_variables <- numeric(0)

while(length(selected_variables) < max_variables){
  if(verbose > 0) print(paste0("On iteration: ", length(selected_variables)+1))
  
  if(verbose > 0) print("Selecting variable")
  cor_vec <- sapply(1:ncol(adt_mat), function(j){
    if(verbose > 1 && ncol(adt_mat) > 10 & j %% floor(ncol(adt_mat)/10) == 0) cat('*')
    tiltedCCA:::.linear_regression(bool_include_intercept = T,
                                   bool_center_x = T,
                                   bool_center_y = T,
                                   bool_scale_x = T,
                                   bool_scale_y = T,
                                   return_type = "r_squared", 
                                   x_mat = reference_dimred,
                                   y_vec = adt_mat[,j])
  })
  candidate_var <- colnames(adt_mat)[which(cor_vec <= cor_threshold)]
  candidate_list[[length(selected_variables)+1]] <- candidate_var
  
  if(verbose > 0) print(paste0(length(candidate_var), " eligble variables"))
  if(length(candidate_var) == 0) break()
  
  idx <- candidate_var[which.max(logpval_vec[candidate_var])]
  selected_variables <- c(selected_variables, idx)
  reference_dimred <- cbind(reference_dimred, adt_mat[,idx])
  colnames(reference_dimred)[ncol(reference_dimred)] <- idx
  
  stopifnot(all(selected_variables %in% colnames(reference_dimred)))
  adt_mat <- adt_mat[,which(!colnames(adt_mat) %in% selected_variables),drop = F]
  if(ncol(adt_mat) == 0) break()
}
panel_3 <- selected_variables

panel_list <- list(most_DE = panel_1,
                   least_aligned = panel_2,
                   no_tcca = panel_3)

save(pbmc, panel_list,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_varSelect_alternatives.RData")

#################
svd_1 <- multiSVD_obj$svd_1
for(i in 1:3){
  print(paste0("Working on custom panel ", i))
  assay_name <- paste0("ADT", i)
  
  adt_mat2 <- pbmc[["ADT"]]@scale.data
  adt_mat2 <- adt_mat2[panel_list[[i]],]
  pbmc[[assay_name]] <- Seurat::CreateAssayObject(counts = adt_mat2)
  pbmc[[assay_name]]@data <- adt_mat2
  pbmc[[assay_name]]@scale.data <- adt_mat2
  pbmc[[assay_name]]@var.features <- rownames(adt_mat2)
  
  adt_mat2 <- t(adt_mat2)
  consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = NULL, mat_2 = adt_mat2,
                                             dims_1 = NULL, dims_2 = 1:ncol(adt_mat2),
                                             dims_consensus = max(ncol(svd_1$u), ncol(adt_mat2)),
                                             svd_1 = svd_1)
  
  set.seed(10)
  umap_res <- Seurat::RunUMAP(consensus_pca$dimred_2)
  umap_mat <- umap_res@cell.embeddings
  rownames(umap_mat) <- colnames(pbmc)
  umap_name <- paste0("adt.umap", i)
  pbmc[[umap_name]] <- Seurat::CreateDimReducObject(umap_mat, assay = "ADT")
  
  set.seed(10)
  umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
  umap_mat <- umap_res@cell.embeddings
  rownames(umap_mat) <- colnames(pbmc)
  umap_name <- paste0("consensusUMAP", i)
  pbmc[[umap_name]] <- Seurat::CreateDimReducObject(umap_mat, assay = "SCT")
}

save(pbmc, panel_list,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_varSelect_alternatives.RData")
