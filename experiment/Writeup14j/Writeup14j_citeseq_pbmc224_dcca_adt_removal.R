rm(list=ls())
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca.RData")
load("../../../../out/Writeup14i/Writeup14i_citeseq_pbmc224_dcca_variable_full_de.RData")
source("../Writeup14f/gene_exploration.R")
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_tmp.RData")
dcca_res <- res
class(dcca_res) <- "dcca"
library(Seurat)

set.seed(10)

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)

summary_mat <- compute_variable_summary(mat = mat_2b, 
                                        common_mat = dcca_decomp$common_mat_2,
                                        metacell_clustering = NA)
for(i in 1:nrow(summary_mat)){
  gene_name <- rownames(summary_mat)[i]
  
  val <- stats::median(sapply(1:length(de_list), function(j){
    idx <- which(rownames(de_list[[j]]) == gene_name)
    if(length(idx) == 0) return(1)
    de_list[[j]][idx, "p_val"]
  }))
  
  summary_mat[i,"kl_div"] <- val
}
summary_mat[,"kl_div"] <- -log10(summary_mat[,"kl_div"])

##################

cutoff_vec <- seq(.7, .5, by = -0.05)

for(cutoff in cutoff_vec){
  print(paste0("Working on cutoff ", cutoff))
  
  set.seed(10)
  pbmc2 <- pbmc
  pbmc2[["adt.umap"]] <- NULL
  pbmc2[["apca"]] <- NULL
  
  # grab which ADT's will be removed
  print("Removing antibodies")
  adt_names <- rownames(summary_mat)[which(summary_mat[,"r_squared"] >= cutoff)]
  
  # compute new ADT PCA
  print("Computing ADT PCA")
  Seurat::DefaultAssay(pbmc2) <- "ADT"
  if(length(adt_names) > 0){
    pbmc2[["ADT"]]@scale.data <- pbmc2[["ADT"]]@scale.data[-which(rownames(pbmc2) %in% adt_names),]
  }
  print(dim(pbmc2[["ADT"]]@scale.data))
  pbmc2 <- Seurat::RunPCA(
    pbmc2, 
    reduction.name = 'apca',
    features = rownames(pbmc2[["ADT"]]@scale.data),
    verbose = F
  )
  
  # compute new ADT umap
  print("Computing ADT UMAP")
  set.seed(10)
  pbmc2 <- Seurat::RunUMAP(
    pbmc2, 
    reduction = 'apca', 
    dims = 1:rank_2, assay = 'ADT', 
    reduction.name = 'adt.umap', 
    reduction.key = 'adtUMAP_'
  )
  
  # compute new multimodal neighbor
  print("Computing WNN")
  pbmc2 <- Seurat::FindMultiModalNeighbors(
    pbmc2, 
    reduction.list = list("pca", "apca"), 
    dims.list = list(1:rank_1, 1:rank_2), 
    modality.weight.name = "RNA.weight"
  )
  
  # compute new WNN UMAP
  print("Computing WNN UMAP")
  set.seed(10)
  pbmc2 <- Seurat::RunUMAP(
    pbmc2, 
    nn.name = "weighted.nn", 
    reduction.name = "wnn.umap", 
    reduction.key = "wnnUMAP_"
  )
  
  # compute consensus PCA's UMAP
  print("Computing Consensus PCA")
  n <- ncol(pbmc2)
  svd_1 <- svd(pbmc2[["pca"]]@cell.embeddings)
  embedding_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
  svd_2 <- svd(pbmc2[["apca"]]@cell.embeddings)
  embedding_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
  embedding_all <- cbind(embedding_1, embedding_2)
  rownames(embedding_all) <- colnames(pbmc2)
  pca_res <- stats::prcomp(embedding_all, center = TRUE, scale. = FALSE)
  consensus_mat <- pca_res$x[,1:min(rank_1, rank_2)] 
  
  set.seed(10)
  consensus_umap <- Seurat::RunUMAP(consensus_mat, 
                                    reduction.key = "umapConsensusPCA_")
  pbmc2[["consensus.umap"]] <- Seurat::CreateDimReducObject(consensus_umap@cell.embeddings,
                                                            assay = "ADT")
  
  # rotate all the UMAPs
  print("Rotating UMAPs")
  anchor_name <- "rna.umap"
  other_names <- c("adt.umap", "wnn.umap", "consensus.umap")
  
  for(umap_name in other_names){
    print(umap_name)
    u_mat1 <- pbmc2[[anchor_name]]@cell.embeddings
    u_mat2 <- pbmc2[[umap_name]]@cell.embeddings
    tmp <- svd(t(u_mat1) %*% u_mat2)
    rotation_mat <- tmp$u %*% t(tmp$v)
    tmp <- u_mat2 %*% t(rotation_mat)
    rownames(tmp) <- rownames(pbmc2@meta.data)
    colnames(tmp) <- colnames(pbmc2[[umap_name]]@cell.embeddings)
    pbmc2[[umap_name]]@cell.embeddings <- tmp
  }
  
  # make all the plots
  print("Plotting")
  reduction_vec <- c(other_names)
  group_vec <- c("celltype.l2")
  main_vec <- c("(ADT)", "(WNN)", "(Consensus PCA)")
  file_vec <- c("adt", "wnn", "consensuspca")
  
  for(i in 1:length(reduction_vec)){
    for(j in 1:length(group_vec)){
      plot1 <- Seurat::DimPlot(pbmc2, reduction = reduction_vec[i],
                               group.by = group_vec[j], label = TRUE,
                               repel = TRUE, label.size = 2.5)
      plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\n", main_vec[i], ", cutoff ", cutoff))
      plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
      ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_adtremoval_cutoff", cutoff, "_", file_vec[i], ".png"),
                      plot1, device = "png", width = 6, height = 5, units = "in")
    }
  }
}