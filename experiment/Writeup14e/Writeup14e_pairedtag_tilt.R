rm(list=ls())
library(Seurat)
library(Signac)

histone_names <- c("H3K4me3", "H3K27me3")

cell_types <- c("Astro_Myoc", "Astro_Nnat", "CA1", "CA23", "CGE",
                "CT", "DG", "Endothelial", "Ependymal", "L23",
                "L4", "L5", "L6", "Microglia", "NP",
                "Oligo_MFOL", "Oligo_MOL", "OPC", "PT", "Pvalb",
                "Sst", "Subiculum")
color_vec <- c("#FF6B2C", "#FF8817", "#5E801F", "#06C65B", "#EA6FE6",
               "#294646", "#00803D", "#B26F13", "#804F0F", "#005D7F",
               "#008FC6", "#05A8EB", "#547180", "#EB8E1A", "#487F80",
               "#EB523A", "#B2402D", "#7F2F23", "#6EC6C6", "#803E7E",
               "#B255B0", "#AAEB32")
color_df <- data.frame(celltype = cell_types, color = color_vec)

cell_list <- list(Non_neurons = c("Astro_Myoc", "Astro_Nnat", "OPC", "Microglia", 
                                  "Oligo_MFOL", "Oligo_MOL", "Endothelial", "Ependymal"),
                  Inhibitory_neurons = c("Pvalb", "Sst", "CGE"),
                  Cortical_neurons = c("CT", "NP", "L6", "L5", "L4", "L23", "PT"),
                  Hippocampal_neurons = c("DG", "Subiculum", "CA1", "CA23"))

#########
weight_val <- c(0.61, 0.69)

for(i in 1:length(histone_names)){
  print(i)
  
  load(paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_",
              histone_names[i], ".RData"))
  
  uniq_celltypes <- sort(unique(pairedtag@meta.data$celltype))
  color_vec <- sapply(uniq_celltypes, function(i){
    color_df[which(color_df$celltype == i),"color"]
  })
  
  plot1 <-  Seurat::VlnPlot(pairedtag, features = "SCT.weight", 
                            group.by = "celltype", sort = TRUE,
                            cols = color_vec,
                            pt.size = 0.1) + Seurat::NoLegend()
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-tag: ", histone_names[i], "\nWNN's RNA weights: Average of ", round(mean(pairedtag$SCT.weight), 2)))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_pairedtag_", histone_names[i], "_wnn_weights.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  #########
  
  embedding_1 <- pairedtag[["pca"]]@cell.embeddings
  embedding_2 <- pairedtag[["lsi"]]@cell.embeddings
  l2_vec_1 <- 1/apply(embedding_1, 1, multiomicCCA:::.l2norm)
  l2_vec_2 <- 1/apply(embedding_2, 1, multiomicCCA:::.l2norm)
  l2_vec_1[is.infinite(l2_vec_1)] <- 0
  l2_vec_2[is.infinite(l2_vec_2)] <- 0
  embedding_1 <- multiomicCCA:::.mult_vec_mat(l2_vec_1, embedding_1)
  embedding_2 <- multiomicCCA:::.mult_vec_mat(l2_vec_2, embedding_2)
  svd_1 <- irlba::irlba(embedding_1, nv = 3)
  svd_2 <- irlba::irlba(embedding_2, nv = 3)
  
  embedding_1b <- embedding_1/svd_1$d[1]*sqrt(nrow(embedding_1))
  embedding_2b <- embedding_2/svd_2$d[1]*sqrt(nrow(embedding_2))
  set.seed(10)
  umap1 <- Seurat::RunUMAP(cbind(embedding_1b, embedding_2b), metric = "euclidean")
  pairedtag[["custom"]] <- Seurat::CreateDimReducObject(umap1@cell.embeddings, key = "custom_", assay = "RNA")
  plot1 <- Seurat::DimPlot(pairedtag, reduction = "custom",
                           group.by = "celltype", label = TRUE,
                           cols = color_vec,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ",  histone_names[i],"\nConsensus-PCA: Weight of 0.5"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_pairedtag_",
                                    histone_names[i],
                                    "_consensus_0.5.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  ###
  
  embedding_1b <- embedding_1/svd_1$d[1]*sqrt(nrow(embedding_1))*weight_val[i]
  embedding_2b <- embedding_2/svd_2$d[1]*sqrt(nrow(embedding_2))*(1-weight_val[i])
  set.seed(10)
  umap2 <- Seurat::RunUMAP(cbind(embedding_1b, embedding_2b), metric = "euclidean")
  pairedtag[["custom"]] <- Seurat::CreateDimReducObject(umap2@cell.embeddings, key = "custom_", assay = "RNA")
  plot1 <- Seurat::DimPlot(pairedtag, reduction = "custom",
                           group.by = "celltype", label = TRUE,
                           cols = color_vec,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ",  histone_names[i],"\nConsensus-PCA: Weight of ", weight_val[i]))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_pairedtag_",
                                    histone_names[i],
                                    "_consensus_", weight_val[i], ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  
  set.seed(10)
  umap3 <- Seurat::RunUMAP(cbind(embedding_1, embedding_2), metric = "cosine")
  pairedtag[["custom"]] <- Seurat::CreateDimReducObject(umap3@cell.embeddings, key = "custom_", assay = "RNA")
  plot1 <- Seurat::DimPlot(pairedtag, reduction = "custom",
                           group.by = "celltype", label = TRUE,
                           cols = color_vec,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ",  histone_names[i],"\nConsensus-PCA: No weighting"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_pairedtag_",
                                    histone_names[i],
                                    "_consensus_false.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
  
  
}


