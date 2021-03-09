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

custom_explained_variance <- function(dcca_decomp, weight_func = function(x){x}){
  stopifnot(class(dcca_decomp) == "dcca_decomp")
  
  n <- nrow(dcca_decomp$common_mat_1)
  p1 <- ncol(dcca_decomp$common_mat_1); p2 <- ncol(dcca_decomp$common_mat_2)
  
  cell_weight_vec1 <- sapply(1:n, function(i){
    val1 <- weight_func(.l2norm(dcca_decomp$common_mat_1[i,]))
    val2 <- weight_func(.l2norm(dcca_decomp$distinct_mat_1[i,]))
    val2/(val1+val2)
  })
  
  cell_weight_vec2 <- sapply(1:n, function(i){
    val1 <- weight_func(.l2norm(dcca_decomp$common_mat_2[i,]))
    val2 <- weight_func(.l2norm(dcca_decomp$distinct_mat_2[i,]))
    val2/(val1+val2)
  })
  
  list(cell_weight_vec1 = cell_weight_vec1,
       cell_weight_vec2 = cell_weight_vec2)
}

custom_explained_variance2 <- function(dcca_decomp, membership_vec){
  stopifnot(class(dcca_decomp) == "dcca_decomp")
  
  n <- nrow(dcca_decomp$common_mat_1)
  mat_1 <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
  weight_vec_1 <- sapply(1:length(unique(membership_vec)), function(i){
    idx <- which(membership_vec == i)
    val <- .l2norm(dcca_decomp$distinct_mat_1[idx,])^2/.l2norm(mat_1[idx,])^2
    max(min(val, 1),0)
  })
  
  mat_2 <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2
  weight_vec_2 <- sapply(1:length(unique(membership_vec)), function(i){
    idx <- which(membership_vec == i)
    val <- .l2norm(dcca_decomp$distinct_mat_2[idx,])^2/.l2norm(mat_2[idx,])^2
    max(min(val, 1),0)
  })
  
  list(weight_vec_1 = weight_vec_1,
       weight_vec_2 = weight_vec_2)
}

pipeline_information_plot <- function(dcca_decomp, true_membership_vec, main = ""){
  tmp <- dcca_information_weight(mat = dcca_decomp$common_mat_1+dcca_decomp$distinct_mat_1,
                                 common_mat = dcca_decomp$common_mat_1)
  weight_vec1 <- tmp$alpha_vec
  tmp <- dcca_information_weight(mat = dcca_decomp$common_mat_2+dcca_decomp$distinct_mat_2,
                                 common_mat = dcca_decomp$common_mat_2)
  weight_vec2 <- tmp$alpha_vec
  weight_mat <- data.frame(Mode_1 = 1-weight_vec1, Mode_2 = 1-weight_vec2, cell_type = true_membership_vec)
  par(mfrow = c(1,1), mar = c(4, 4, 3.5, 0.5))
  information_plot(weight_mat, col_vec = 1:length(unique(true_membership_vec)), reorder_types = F,
                   main = main)
  title(xlab="Cell-type", mgp=c(1,1,0))
}

