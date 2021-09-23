compute_metacells <- function(seurat_obj,
                              embedding_1,
                              embedding_2,
                              target_num_clusters = round(nrow(embedding_1)/100),
                              starting_resolution = 1,
                              scaling_factor = 1.5,
                              verbose = T){
  stopifnot(nrow(embedding_1) == nrow(embedding_2))
  
  all_mat <- cbind(embedding_1, embedding_2)
  snn <- Seurat::FindNeighbors(all_mat,
                               compute.SNN = T)
  
  snn$snn@assay.used <- "RNA"
  seurat_obj[["newgraph"]] <- snn$snn
  # (length(pbmc@graphs$snn@x)-pbmc@graphs$snn@Dim[1])/2
  
  resolution <- starting_resolution
  while(TRUE){
    seurat_obj2 <- Seurat::FindClusters(seurat_obj, 
                                        graph.name = "newgraph",
                                        resolution = resolution)
    num_clusters <- length(unique(seurat_obj2@meta.data[,"seurat_clusters"]))
    if(num_clusters >= target_num_clusters) break()
    
    if(verbose) print(paste0("Trying resolution ", resolution, 
                             ": Number of clusters ", num_clusters))
    resolution <- round(resolution * scaling_factor)
  }
 
  as.factor(seurat_obj2@meta.data[,"seurat_clusters"])
}