rm(list=ls())

histone_names <- c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3")

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


for(i in 1:length(histone_names)){
  print(i)
  
  load(paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_",
              histone_names[i], ".RData"))
  load(paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_",
              histone_names[i], ".RData"))
  
  # "binarizing"
  dat_1b <- dat_1
  dat_1b@x <- rep(1, length(dat_1b@x))
  
  # Cells with less than 200 features in both DNA and RNA matrices were removed. 
  keep_vec <- rep(1, ncol(dat_1b))
  colsum_vec <- sparseMatrixStats::colSums2(dat_1b)
  quantile(colsum_vec)
  keep_vec[colsum_vec < 200] <- 0
  keep_vec[which(pairedtag@meta.data$nFeature_RNA <= 200)] <- 0
  table(keep_vec)
  pairedtag[["keep"]] <- keep_vec
  
  # "removing the 5% highest covered bins"
  rowsum_vec <- sparseMatrixStats::rowSums2(dat_1b[,which(keep_vec == 1)])
  idx1 <- which(rowsum_vec >= quantile(rowsum_vec, probs = 0.95))
  idx2 <- which(rowsum_vec <= 4)
  dat_1b <- dat_1b[-c(unique(c(idx1, idx2))),]
  dat_1 <- dat_1[-c(unique(c(idx1, idx2))),]
  
  pairedtag[["DNA"]] <- Seurat::CreateAssayObject(counts = dat_1)
  Seurat::DefaultAssay(pairedtag) <- "DNA"
  pairedtag[["DNA"]]@data <- pairedtag[["DNA"]]@counts
  pairedtag <- subset(pairedtag, keep == 1)
  set.seed(10)
  dat_1c <- Matrix::t(pairedtag[["DNA"]]@data)
  rowsum_vec <- sparseMatrixStats::rowSums2(dat_1c)
  diag_mat <- Matrix::Diagonal(x = 1/rowsum_vec*100)
  dat_1c <- diag_mat %*% dat_1c
  center_vec <- sparseMatrixStats::colMeans2(dat_1c)
  sd_vec <- sparseMatrixStats::colSds(dat_1c)
  svd_res <- irlba::irlba(dat_1c, nv = 50, 
                          scale = sd_vec, 
                          center = center_vec)
  dim_red <- multiomicCCA:::.mult_mat_vec(svd_res$u, svd_res$d)
  dim_red <- scale(dim_red, center = T, scale = T)
  rownames(dim_red) <- rownames(pairedtag@meta.data)
  colnames(dim_red) <- paste0("lsi_", 1:50)
  pairedtag[["lsi"]] <- Seurat::CreateDimReducObject(embedding = dim_red, 
                                                     key = "lsi_", assay = "DNA")
  
  set.seed(10)
  pairedtag <- Seurat::RunUMAP(pairedtag, 
                               reduction="lsi", 
                               dims=1:25,
                               metric = "euclidean",
                               reduction.name="umap.dna", 
                               reduction.key="dnaUMAP_")
  
  uniq_celltypes <- sort(unique(pairedtag@meta.data$celltype))
  color_vec <- sapply(uniq_celltypes, function(i){
    color_df[which(color_df$celltype == i),"color"]
  })
  
  plot1 <- Seurat::DimPlot(pairedtag, reduction = "umap.dna",
                           group.by = "celltype", label = TRUE,
                           cols = color_vec,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("Paired-Tag: ",  histone_names[i]," (Histone)"))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_",
                                    histone_names[i],
                                    "_histone.png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}

# 
# png("../../../../out/figures/Writeup14d/Writeup14d_pairedtag_H3K27me3_histone_heatmap.png",
#     height = 3000, width = 3000, units = "px", res = 300)
# plot_scores_heatmap(dim_red, membership_vec = as.factor(pairedtag@meta.data$celltype),
#                     bool_log = T, bool_center = F, bool_scale = F,
#                     bool_rownormalize_before = F,
#                     bool_rownormalize_after = F,
#                     scaling_power = 5)
# graphics.off()
