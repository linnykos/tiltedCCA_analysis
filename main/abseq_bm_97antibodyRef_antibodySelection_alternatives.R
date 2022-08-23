rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/abseq_bm97Ref_differential.RData")
load("../../../out/main/abseq_bm97Ref_tcca.RData")
source("bm_97antibodyRef_colorPalette.R")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# gather all the relevant ingredients
Seurat::DefaultAssay(bm) <- "AB"
antibody_names <- rownames(bm)
logpval_vec <- sapply(1:length(antibody_names), function(k){
  if(k %% floor(length(antibody_names)/10) == 0) cat('*')
  antibody <- antibody_names[k]
  
  # cycle through all the celltypes
  celltype_vec <- sapply(1:length(adt_de_list$level_vec), function(i){
    idx <- which(adt_de_list$combn_mat == i, arr.ind = T)[,2]
    vec <-  sapply(idx, function(j){
      idx <- which(rownames(adt_de_list$de_list[[j]]) == antibody)
      if(length(idx) == 0) return(1)
      adt_de_list$de_list[[j]][idx, "p_val"]
    })
    stats::quantile(vec, probs = 0.75)
  })
  
  max(-log10(celltype_vec))
})
logpval_vec <- pmin(logpval_vec, 300)
names(logpval_vec) <- rownames(bm)

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(multiSVD_obj, 
                                                    bool_modality_1_full = F,
                                                    bool_modality_2_full = F,
                                                    verbose = 0)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
reference_dimred <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_dimred")
adt_mat <- t(as.matrix(bm[["AB"]]@scale.data))

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
num_variables <- 10

# First panel: Select the antibodies with the highest p-vales
panel_1 <- names(logpval_vec)[order(logpval_vec, decreasing = T)[1:num_variables]]

# Second panel: Select the antibodies deviating the most from the common space
panel_2 <- names(rsquare_vec)[order(rsquare_vec, decreasing = F)[1:num_variables]]

# Third panel: Select the antibodies using the proposed strategy, but only consider the deviation from the RNA subspace
reference_dimred <- tiltedCCA:::.mult_mat_vec(multiSVD_obj$svd_1$u, multiSVD_obj$svd_1$d)
adt_mat <- t(as.matrix(bm[["AB"]]@scale.data))
logpval_vec <- logpval_vec[colnames(adt_mat)]

res <- tiltedCCA:::.generic_variable_selection(
  bool_maximizing = T, 
  cor_threshold = 0.85,
  initial_mat = reference_dimred, # can be NULL
  mat = adt_mat,
  num_variables = num_variables,
  return_candidate_list = F,
  vec = logpval_vec,
  verbose = 2
)
panel_3 <- res

panel_list <- list(most_DE = panel_1,
                   least_aligned = panel_2,
                   no_tcca = panel_3)

###########################

adt_mat <- t(as.matrix(bm[["AB"]]@scale.data))
correlation_vec <- sapply(1:ncol(adt_mat), function(j){
  print(j)
   
  df <- data.frame(cbind(adt_mat[,j], multiSVD_obj$svd_1$u))
  colnames(df)[1] <- "y"
  lm_res <- stats::lm(y ~ ., data = df)
  summary(lm_res)$r.squared
})
names(correlation_vec) <- colnames(adt_mat)
correlation_vec <- correlation_vec[!is.nan(correlation_vec)]
panel_list[["most_uncorrelated"]] <- names(correlation_vec)[order(correlation_vec, decreasing = F)[1:10]]

save(bm, panel_list,
     date_of_run, session_info,
     file = "../../../out/main/abseq_bm97Ref_varSelect_alternatives.RData")

#################
svd_1 <- multiSVD_obj$svd_1
for(i in 1:4){
  print(paste0("Working on custom panel ", i))
  assay_name <- paste0("AB", i)
  
  adt_mat2 <- bm[["AB"]]@scale.data
  adt_mat2 <- adt_mat2[panel_list[[i]],]
  bm[[assay_name]] <- Seurat::CreateAssayObject(counts = adt_mat2)
  bm[[assay_name]]@data <- adt_mat2
  bm[[assay_name]]@scale.data <- adt_mat2
  bm[[assay_name]]@var.features <- rownames(adt_mat2)
  
  adt_mat2 <- t(adt_mat2)
  set.seed(10)
  consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = NULL, mat_2 = adt_mat2,
                                             dims_1 = NULL, dims_2 = 1:ncol(adt_mat2),
                                             dims_consensus = 1:20,
                                             svd_1 = svd_1)
  dimred_name <- paste0("consensusPCA", i)
  bm[[dimred_name]] <- Seurat::CreateDimReducObject(consensus_pca$dimred_consensus, 
                                                    assay = "RNA",
                                                    key = paste0("cPC", i))
  
  set.seed(10)
  umap_res <- Seurat::RunUMAP(consensus_pca$dimred_2)
  umap_mat <- umap_res@cell.embeddings
  rownames(umap_mat) <- colnames(bm)
  umap_name <- paste0("adt.umap", i)
  bm[[umap_name]] <- Seurat::CreateDimReducObject(umap_mat, assay = "RNA")
  
  set.seed(10)
  umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
  umap_mat <- umap_res@cell.embeddings
  rownames(umap_mat) <- colnames(bm)
  umap_name <- paste0("consensusUMAP", i)
  bm[[umap_name]] <- Seurat::CreateDimReducObject(umap_mat, assay = "RNA")
}

save(bm, panel_list, 
     date_of_run, session_info,
     file = "../../../out/main/abseq_bm97Ref_varSelect_alternatives.RData")
