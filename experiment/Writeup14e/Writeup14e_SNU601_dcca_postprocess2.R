rm(list=ls())
load("../../../../out/Writeup14e/Writeup14e_SNU_dcca.RData")

cna_embeddings <- multiomicCCA:::.prepare_embeddings(dcca_res, data_1 = T, data_2 = F, 
                                 center = F, 
                                 renormalize = F)
atac_embeddings <- multiomicCCA:::.prepare_embeddings(dcca_res, data_1 = F, data_2 = T, 
                                                     center = F, 
                                                     renormalize = F)
atac_embeddings <- multiomicCCA:::.normalize_embeddings(atac_embeddings, normalization_type = "signac_itself")

for(i in 1:3){
  l2_vec <- 1/apply(cna_embeddings[[i]], 1, multiomicCCA:::.l2norm)
  l2_vec[is.infinite(l2_vec)] <- 0
  cna_embeddings[[i]] <- multiomicCCA:::.mult_vec_mat(l2_vec,  cna_embeddings[[i]])
  
  l2_vec <- 1/apply(atac_embeddings[[i]], 1, multiomicCCA:::.l2norm)
  l2_vec[is.infinite(l2_vec)] <- 0
  atac_embeddings[[i]] <- multiomicCCA:::.mult_vec_mat(l2_vec,  atac_embeddings[[i]])
}

combined_common_embedding <- cbind(cna_embeddings[["common"]], atac_embeddings[["common"]])

set.seed(10)
tmp <- Seurat::RunUMAP(combined_common_embedding)@cell.embeddings
rownames(tmp) <- colnames(SNU[["CNA"]]@counts)
SNU[["combined"]] <- Seurat::CreateDimReducObject(embeddings = tmp,
                                                        key = "consensusUMAP_", assay = "ATAC")

set.seed(10)
tmp <- Seurat::RunUMAP(cna_embeddings[["distinct"]])@cell.embeddings
rownames(tmp) <- colnames(SNU[["CNA"]]@counts)
SNU[["distinct"]] <- Seurat::CreateDimReducObject(embeddings = tmp,
                                                  key = "distinct_", assay = "CNA")

set.seed(10)
tmp <- Seurat::RunUMAP(atac_embeddings[["distinct"]])@cell.embeddings
rownames(tmp) <- colnames(SNU[["CNA"]]@counts)
SNU[["distinct2"]] <- Seurat::CreateDimReducObject(embeddings = tmp,
                                                  key = "distinct2_", assay = "ATAC")

set.seed(10)
tmp <- Seurat::RunUMAP(cna_embeddings[["everything"]])@cell.embeddings
rownames(tmp) <- colnames(SNU[["CNA"]]@counts)
SNU[["everything"]] <- Seurat::CreateDimReducObject(embeddings = tmp,
                                                  key = "distinct_", assay = "CNA")

set.seed(10)
tmp <- Seurat::RunUMAP(atac_embeddings[["everything"]])@cell.embeddings
rownames(tmp) <- colnames(SNU[["CNA"]]@counts)
SNU[["everything2"]] <- Seurat::CreateDimReducObject(embeddings = tmp,
                                                   key = "distinct2_", assay = "ATAC")


##############

group_vec <- c("clone") #, "metacell")

for(j in 1:length(group_vec)){
  plot1 <- Seurat::DimPlot(SNU, reduction = "distinct",
                           group.by = group_vec[j], label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (CNA)\nDistinct view"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_dcca_cna_distinct_umap_custom_", group_vec[j], ".png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
  
  plot1 <- Seurat::DimPlot(SNU, reduction = "distinct2",
                           group.by = group_vec[j], label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (ATAC)\nDistinct view"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_dcca_atac_distinct_umap_custom_", group_vec[j], ".png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
  plot1 <- Seurat::DimPlot(SNU, reduction = "everything",
                           group.by = group_vec[j], label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (CNA)\nEverything view"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_dcca_cna_everything_umap_custom_", group_vec[j], ".png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
  
  plot1 <- Seurat::DimPlot(SNU, reduction = "everything2",
                           group.by = group_vec[j], label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (ATAC)\nEverything view"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_dcca_atac_everything_umap_custom_", group_vec[j], ".png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
  
  plot1 <- Seurat::DimPlot(SNU, reduction = "combined",
                           group.by = group_vec[j], label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU\nBoth Common"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_dcca_bothcommon_umap_custom_", group_vec[j], ".png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
}
