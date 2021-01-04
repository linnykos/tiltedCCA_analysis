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

# apply_seurat_map <- function(mat, seurat_obj){
#   rownames(mat) <- colnames(pbmc)
#   umap_res <- Seurat::RunUMAP(mat, reduction.key = 'tmp')
#   umap_res@cell.embeddings
# }