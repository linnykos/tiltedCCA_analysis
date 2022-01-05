rm(list=ls())
library(Seurat)
load("../../../../out/Writeup14h/Writeup14h_citeseq_pbmc224_dcca5_enumeration.RData")
discretization_gridsize <- rev(seq(0, 1, length = length(dcca_list)))

for(kk in 1:length(discretization_gridsize)){
  tilt <- discretization_gridsize[kk]
  pbmc2 <- pbmc
  print(paste0("Tilt: ", tilt))
  
  # common
  dcca_res <- dcca_list[[kk]]
  dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)
  
  n <- nrow(pbmc2@meta.data)
  svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_1, 
                                         K = rank_1, K_full_rank = F,
                                         rescale = F,
                                         mean_vec = F, sd_vec = F,
                                         symmetric = F)
  svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$common_mat_2, K = rank_2, K_full_rank = F,
                                         rescale = F,
                                         mean_vec = F, sd_vec = F,
                                         symmetric = F)
  dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d/svd_1$d[1]*sqrt(n))
  l2_vec <- apply(dimred_1, 1, function(x){multiomicCCA:::.l2norm(x)})
  dimred_1 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_1)
  dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d/svd_2$d[1]*sqrt(n))
  l2_vec <- apply(dimred_2, 1, function(x){multiomicCCA:::.l2norm(x)})
  dimred_2 <- multiomicCCA:::.mult_vec_mat(1/l2_vec, dimred_2)
  
  set.seed(10)
  common_umap <- Seurat::RunUMAP(cbind(dimred_1, dimred_2), 
                                 metric = "euclidean",
                                 reduction.key = "umapCommon_")
  rownames(common_umap@cell.embeddings) <- rownames(pbmc2@meta.data)
  pbmc2[["dcca_common"]] <- Seurat::CreateDimReducObject(common_umap@cell.embeddings, 
                                                         assay = "SCT")
  
  
  #######################################
  
  # distinct 1
  svd_1 <- multiomicCCA:::.svd_truncated(dcca_decomp$distinct_mat_1, 
                                         K = rank_1, K_full_rank = F,
                                         rescale = F,
                                         mean_vec = F, sd_vec = F,
                                         symmetric = F)
  dimred_1 <- multiomicCCA:::.mult_mat_vec(svd_1$u, svd_1$d)
  
  set.seed(10)
  distinct1_umap <- Seurat::RunUMAP(dimred_1, 
                                    metric = "cosine",
                                    reduction.key = "umapDistinct1_")
  rownames(distinct1_umap@cell.embeddings) <- rownames(pbmc2@meta.data)
  pbmc2[["dcca_distinct1"]] <- Seurat::CreateDimReducObject(distinct1_umap@cell.embeddings, 
                                                            assay = "SCT")
  
  #######################################
  
  # distinct 2
  svd_2 <- multiomicCCA:::.svd_truncated(dcca_decomp$distinct_mat_2, 
                                         K = rank_2, K_full_rank = F,
                                         rescale = F,
                                         mean_vec = F, sd_vec = F,
                                         symmetric = F)
  dimred_2 <- multiomicCCA:::.mult_mat_vec(svd_2$u, svd_2$d)
  
  set.seed(10)
  distinct2_umap <- Seurat::RunUMAP(dimred_2, 
                                    metric = "cosine",
                                    reduction.key = "umapDistinct2_")
  rownames(distinct2_umap@cell.embeddings) <- rownames(pbmc2@meta.data)
  pbmc2[["dcca_distinct2"]] <- Seurat::CreateDimReducObject(distinct2_umap@cell.embeddings,
                                                            assay = "ADT")
  
  #######################################
  
  # rotate
  anchor_name <- "rna.umap"
  other_names <- c("dcca_common", "dcca_distinct1", "dcca_distinct2")
  
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
  
  # plot according to clones
  reduction_vec <- other_names
  group_vec <- c("celltype.l2")
  main_vec <- c(paste0("(D-CCA, Common, Tilt: ", tilt, ")"), 
                paste0("(D-CCA, Distinct 1, Tilt: ", tilt, ")"), 
                paste0("(D-CCA, Distinct 2, Tilt: ", tilt, ")"))
  file_vec <- c("dcca-common", "dcca-distinct1", "dcca-distinct2")
  
  for(i in 1:length(reduction_vec)){
    for(j in 1:length(group_vec)){
      plot1 <- Seurat::DimPlot(pbmc2, reduction = reduction_vec[i],
                               group.by = group_vec[j], label = TRUE,
                               repel = TRUE, label.size = 2.5)
      plot1 <- plot1 + ggplot2::ggtitle(paste0("Human PBMC (Cite-seq):\n", main_vec[i]))
      plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
      ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14h/Writeup14h_citeseq_pbmc224_dcca5_tilt", tilt, "_", file_vec[i], ".png"),
                      plot1, device = "png", width = 6, height = 5, units = "in")
    }
  }
}
