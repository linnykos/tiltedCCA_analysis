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

plot_scores <- function(dat, mode = 1){
  if(mode == 1){
    zlim <- range(c(dat$common_score, dat$distinct_score_1, dat$distinct_score_2))
    par(mfrow = c(1,3))
    image(t(dat$common_score), zlim = zlim, main = "Common score")
    image(t(dat$distinct_score_1), zlim = zlim, main = paste0("Distinct score 1"))
    image(t(dat$distinct_score_2), zlim = zlim, main = paste0("Distinct score 2"))
    
  } else {
    zlim <- range(c(dat$common_score_1, dat$common_score_2, dat$distinct_score_1, dat$distinct_score_2))
    par(mfrow = c(1,4))
    image(t(dat$common_score_1), zlim = zlim, main = "Common score 1")
    image(t(dat$common_score_2), zlim = zlim, main = "Common score 2")
    image(t(dat$distinct_score_1), zlim = zlim, main = paste0("Distinct score 1"))
    image(t(dat$distinct_score_2), zlim = zlim, main = paste0("Distinct score 2"))
    
  }
  
  invisible()
}

plot_embeddings <- function(dat, membership_vec){
  par(mfrow = c(1,3))
  set.seed(10)
  tmp <- extract_embedding(dat, common_1 = T, common_2 = T, distinct_1 = F, distinct_2 = F, mode = "dmca")
  plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Common view",
       xlab = "UMAP 1", ylab = "UMAP 2")
  set.seed(10)
  tmp <- extract_embedding(dat, common_1 = F, common_2 = F, distinct_1 = T, distinct_2 = T, mode = "dmca")
  plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Distinct views",
       xlab = "UMAP 1", ylab = "UMAP 2")
  set.seed(10)
  tmp <- extract_embedding(dat, common_1 = T, common_2 = T, distinct_1 = T, distinct_2 = T, mode = "dmca")
  plot(tmp[,1], tmp[,2], asp = T, pch = 16, col = membership_vec, main = "Entire view",
       xlab = "UMAP 1", ylab = "UMAP 2")
  
  invisible()
}


