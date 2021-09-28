wnn_custom <- function(mat_1,
                       mat_2,
                       svd_1,
                       svd_2,
                       dims_1,
                       dims_2){
  seurat_obj <- Seurat::CreateSeuratObject(counts = t(mat_1))
  seurat_obj[["ATAC"]] <- Seurat::CreateAssayObject(counts = t(mat_2))
  
  dimred1 <- multiomicCCA:::.mult_mat_vec(svd_1$u[,dims_1,drop = F], svd_1$d[dims_1])
  colnames(dimred1) <- paste0("dimred1_", dims_1)
  seurat_obj[["dimred1"]] <- Seurat::CreateDimReducObject(embedding = dimred1, 
                                                          key = "dimred1_", assay = "RNA")
  
  dimred2 <- multiomicCCA:::.mult_mat_vec(svd_2$u[,dims_2,drop = F], svd_2$d[dims_2])
  colnames(dimred2) <- paste0("dimred2_", dims_2)
  seurat_obj[["dimred2"]] <- Seurat::CreateDimReducObject(embedding = dimred2, 
                                                          key = "dimred2_", assay = "ATAC")
  
  seurat_obj <- Seurat::FindMultiModalNeighbors(
    seurat_obj, 
    reduction.list = list("dimred1", "dimred2"), 
    dims.list = list(1:length(dims_1), 1:length(dims_2))
  )
  
  seurat_obj <- Seurat::RunUMAP(seurat_obj, 
                                nn.name = "weighted.nn", 
                                reduction.name = "wnn.umap", 
                                reduction.key = "wnnUMAP_")
  
  seurat_obj[["wnn.umap"]]@cell.embeddings
}