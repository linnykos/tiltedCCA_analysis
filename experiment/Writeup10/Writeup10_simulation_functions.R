analyze_seurat_pipeline <- function(mat_1, mat_2){
  seurat_obj <- form_seurat_obj(mat_1, mat_2)
  Seurat::DefaultAssay(seurat_obj) <- "mode1"
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = rownames(seurat_obj), 
                                  do.scale = F, verbose = F)
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = rownames(seurat_obj), 
                               npcs = K, verbose = F, reduction.key = "PC1_", reduction.name = "pca1")
  seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = 'pca1', dims = 1:K, assay = 'mode1', 
                                reduction.name = 'umap1', reduction.key = 'UMAP1_', verbose = F)
  Seurat::DefaultAssay(seurat_obj) <- "mode2"
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = rownames(seurat_obj), 
                                  do.scale = F, verbose = F)
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = rownames(seurat_obj), 
                               npcs = K, verbose = F, reduction.key = "PC2_", reduction.name = "pca2") 
  seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = 'pca2', dims = 1:K, assay = 'mode2', 
                                reduction.name = 'umap2', reduction.key = 'UMAP2_', verbose = F)
  seurat_obj <- Seurat::FindMultiModalNeighbors(seurat_obj, reduction.list = list("pca1", "pca2"), 
                                                weighted.nn.name = "weighted.nn", dims.list = list(1:K, 1:K), 
                                                modality.weight.name = "mode1.weight", verbose = F)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                                reduction.key = "wnnUMAP_", verbose = F)
  
  seurat_obj
}

construct_noiseless_data <- function(common_score, distinct_score_1, distinct_score_2,
                                     coef_mat_1, coef_mat_2){
  
  rank_c <- ncol(common_score)
  common_mat_1 <- common_score %*% coef_mat_1[1:rank_c,,drop = F]
  distinct_mat_1 <- distinct_score_1 %*% coef_mat_1
  common_mat_2 <- common_score %*% coef_mat_2[1:rank_c,,drop = F] 
  distinct_mat_2 <- distinct_score_2 %*% coef_mat_2
  
  structure(list(common_score = common_score,
      distinct_score_1 = distinct_score_1,
      distinct_score_2 = distinct_score_2,
      common_mat_1 = common_mat_1, common_mat_2 = common_mat_2, 
      distinct_mat_1 = distinct_mat_1, distinct_mat_2 = distinct_mat_2), class = "dcca_decomp")
}