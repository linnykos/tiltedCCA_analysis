rm(list=ls())
load("../../../../out/Writeup14e/Writeup14e_SNU_dcca.RData")
SNU[["metacell"]] <- dcca_res$metacell_clustering

set.seed(10)
tmp <- Seurat::RunUMAP(dcca_res$consensus_pca_mat)@cell.embeddings
rownames(tmp) <- colnames(SNU[["CNA"]]@counts)
SNU[["umap.consensus"]] <- Seurat::CreateDimReducObject(embeddings = tmp,
                                                        key = "consensusUMAP_", assay = "ATAC")


# first plot according to clones
reduction_vec <- c("umap.atac", "umap.cna", "wnn.umap", "umap.consensus")
group_vec <- c("clone", "metacell")
main_vec <- c("(ATAC)", "(Copy number)", "(WNN)", "(Consensus PCA)")
file_vec <- c("atac", "cna", "wnn", "consensus")
  
for(i in 1:length(reduction_vec)){
  for(j in 1:length(group_vec)){
    plot1 <- Seurat::DimPlot(SNU, reduction = reduction_vec[i],
                             group.by = group_vec[j], label = TRUE,
                             repel = TRUE, label.size = 2.5)
    plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU ", main_vec[i], ": ", group_vec[j]))
    plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
    ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_", file_vec[i], "_", group_vec[j], ".png"),
                    plot1, device = "png", width = 5, height = 5, units = "in")
  }
}

###############3

png(paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_dcca_scores.png"),
    height = 1200, width = 2500, units = "px", res = 300)
par(mfrow = c(1,3), mar = c(4,4,4,0.5))
multiomicCCA::plot_scores_heatmap.dcca(dcca_res,
                                       membership_vec = as.factor(SNU@meta.data$clone),
                                       log_scale = T, scaling_power = 2)
graphics.off()

#######################

SNU[["common"]] <- Seurat::CreateDimReducObject(embedding = cna_embeddings[[1]], key = "common_", assay = "CNA")
SNU[["distinct"]] <- Seurat::CreateDimReducObject(embedding = cna_embeddings[[2]], key = "distinct_", assay = "CNA")
SNU[["everything"]] <- Seurat::CreateDimReducObject(embedding = cna_embeddings[[3]], key = "everything_", assay = "CNA")

SNU[["common2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[1]], key = "common2_", assay = "ATAC")
SNU[["distinct2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[2]], key = "distinct2_", assay = "ATAC")
SNU[["everything2"]] <- Seurat::CreateDimReducObject(embedding = atac_embeddings[[3]], key = "everything2_", assay = "ATAC")

SNU[["combined"]] <- Seurat::CreateDimReducObject(embedding = combined_common_umap, key = "combined_", assay = "CNA")

title_vec <- c("Common view", "Distinct view", "Everything view")
main_vec <- c("common", "distinct", "everything")

for(i in 1:3){
  plot1 <- Seurat::DimPlot(SNU, reduction = main_vec[i],
                           group.by = "clone", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (CNA)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_dcca_cna_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
}

# plot ATAC embeddings
for(i in 1:3){
  plot1 <- Seurat::DimPlot(SNU, reduction = paste0(main_vec[i], "2"),
                           group.by = "clone", label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (ATAC)\n", title_vec[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_dcca_atac_", main_vec[i], "_umap.png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
}

# plot the combined view
plot1 <- Seurat::DimPlot(SNU, reduction = "combined",
                         group.by = "clone", label = TRUE,
                         repel = TRUE, label.size = 2.5)
plot1 <- plot1 + ggplot2::ggtitle("SNU\nBoth Common")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_dcca_bothcommon_umap.png"),
                plot1, device = "png", width = 5, height = 5, units = "in")

