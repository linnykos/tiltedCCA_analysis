rm(list=ls())
load("../../../out/main/citeseq_pbmc224_differential.RData")
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")
source("pbmc_224antibody_colorPalette.R")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

Seurat::DefaultAssay(pbmc) <- "ADT"
logpval_vec <- sapply(rownames(pbmc), function(gene){
  val <- stats::quantile(
    sapply(1:length(adt_de_list), function(j){
      idx <- which(rownames(adt_de_list[[j]]) == gene)
      if(length(idx) == 0) return(1)
      adt_de_list[[j]][idx, "p_val"]
    }), probs = 0.75
  )
  
  -log10(val)
})
# tmp <- sapply(1:length(adt_de_list), function(j){
#   adt_de_list[[j]][which(rownames(adt_de_list[[j]]) == "CD52"), "p_val"]
# })
names(logpval_vec) <- rownames(pbmc)
quantile(logpval_vec)
max_val <- max(logpval_vec[which(!is.infinite(logpval_vec))])
logpval_vec <- pmin(logpval_vec, 300)

set.seed(10)
variable_selection_res <- tiltedCCA:::postprocess_variable_selection(
  input_obj = multiSVD_obj,
  input_mat = t(as.matrix(pbmc[["ADT"]]@scale.data)),
  logpval_vec = logpval_vec,
  cor_threshold = 0.85,
  input_assay = 2,
  max_variables = 5,
  min_subsample_cell = 5000,
  seurat_celltype_variable = "celltype.l2",
  seurat_obj = pbmc,
  verbose = 2
)

# quantile(logpval_vec)
# logpval_vec[variable_selection_res$selected_variables]
# distinct_mat <- multiSVD_obj$distinct_mat_2
# cor_mat <- stats::cor(distinct_mat[,variable_selection_res$selected_variables])

adt_mat2 <- pbmc[["ADT"]]@scale.data
adt_mat2 <- adt_mat2[variable_selection_res$selected_variables,]
pbmc[["ADT2"]] <- Seurat::CreateAssayObject(counts = adt_mat2)
pbmc[["ADT2"]]@data <- adt_mat2
pbmc[["ADT2"]]@scale.data <- adt_mat2
pbmc[["ADT2"]]@var.features <- rownames(adt_mat2)

Seurat::DefaultAssay(pbmc) <- "ADT2"
pbmc <- Seurat::RunPCA(pbmc, reduction.name = 'apca2', 
                       npcs = nrow(adt_mat2), 
                       verbose = F)

set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, reduction = 'apca2', dims = 1:(ncol(pbmc[["apca2"]]@cell.embeddings)-1), 
                        assay = 'ADT2',
                        reduction.name = 'adt2.umap', reduction.key = 'adt2UMAP_')

set.seed(10)
pbmc <- Seurat::FindMultiModalNeighbors(pbmc, 
                                        reduction.list = list("pca", "apca2"),
                                        dims.list = list(1:40, 1:(ncol(pbmc[["apca2"]]@cell.embeddings)-1)), 
                                        modality.weight.name = "RNA.weight",
                                        weighted.nn.name = "weighted.nn2")
set.seed(10)
pbmc <- Seurat::RunUMAP(pbmc, nn.name = "weighted.nn2", 
                        reduction.name = "wnn2.umap", 
                        reduction.key = "wnn2UMAP_")

svd_1 <- multiSVD_obj$svd_1
adt_mat2 <- t(adt_mat2)
consensus_pca <- tiltedCCA:::consensus_pca(mat_1 = NULL, mat_2 = adt_mat2,
                                           dims_1 = NULL, dims_2 = 1:ncol(adt_mat2),
                                           dims_consensus = 1:max(ncol(svd_1$u), ncol(adt_mat2)),
                                           svd_1 = svd_1, verbose = 1)
save(variable_selection_res, pbmc, 
     logpval_vec, 
     consensus_pca,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_varSelect.RData")

pbmc[["consensusPCA"]] <- Seurat::CreateDimReducObject(consensus_pca$dimred_consensus, 
                                                       assay = "SCT",
                                                       key = "cPC")

set.seed(10)
umap_res <- Seurat::RunUMAP(consensus_pca$dimred_consensus)
umap_mat <- umap_res@cell.embeddings
rownames(umap_mat) <- colnames(pbmc)
pbmc[["consensusUMAP"]] <- Seurat::CreateDimReducObject(umap_mat, 
                                                        assay = "SCT",
                                                        key = "cUMAP")

save(variable_selection_res, pbmc, 
     logpval_vec, 
     consensus_pca,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_varSelect.RData")

######################3

plot2 <- Seurat::DimPlot(pbmc, reduction = "adt2.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot2 <- plot2 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nADT: Selected antibodies, Thres: ", variable_selection_res$cor_threshold))
plot2 <- plot2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_varSelect_adt-umap.png"),
                plot2, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(pbmc, reduction = "wnn2.umap",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nWNN: Selected antibodies, Thres: ", variable_selection_res$cor_threshold))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_varSelect_wnn-umap.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

plot3 <- Seurat::DimPlot(pbmc, reduction = "consensusUMAP",
                         group.by = "celltype.l2", label = TRUE,
                         repel = TRUE, label.size = 2.5,
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + ggplot2::ggtitle(paste0("PBMC (CITE-Seq, RNA+224 ADT)\nConsensusPCA: Selected antibodies, Thres: ", variable_selection_res$cor_threshold))
plot3 <- plot3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_varSelect_consensusPCA-umap.png"),
                plot3, device = "png", width = 6, height = 5, units = "in")

####

plot2 <- Seurat::DimPlot(pbmc, reduction = "adt2.umap",
                         group.by = "celltype.l2",
                         cols = col_palette,
                         raster = FALSE)
plot2 <- plot2 + Seurat::NoLegend() + Seurat::NoAxes()
plot2 <- plot2 + ggplot2::ggtitle("")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_varSelect_adt-umap_cleaned.png"),
                plot2, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot3 <- Seurat::DimPlot(pbmc, reduction = "wnn2.umap",
                         group.by = "celltype.l2",
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + Seurat::NoLegend() + Seurat::NoAxes()
plot3 <- plot3 + ggplot2::ggtitle("")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_varSelect_wnn-umap_cleaned.png"),
                plot3, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

plot3 <- Seurat::DimPlot(pbmc, reduction = "consensusUMAP",
                         group.by = "celltype.l2",
                         cols = col_palette,
                         raster = FALSE)
plot3 <- plot3 + Seurat::NoLegend() + Seurat::NoAxes()
plot3 <- plot3 + ggplot2::ggtitle("")
ggplot2::ggsave(filename = paste0("../../../out/figures/main/citeseq_pbmc224_varSelect_consensusPCA-umap_cleaned.png"),
                plot3, device = "png", width = 3, height = 3, units = "in",
                dpi = 500)

################################

rsquare_vec <- variable_selection_res$cor_vec_intial
rsquare_vec <- rsquare_vec[names(logpval_vec)]
png("../../../out/figures/main/citeseq_pbmc224_differential_protein.png",
    height = 3500, width = 2500, res = 500, units = "px")
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "PBMC (CITE-Seq, RNA+224 ADT)\nAntibody differentiability vs. alignment",
                           bool_mark_ymedian = F,
                           bool_polygon_mean = F,
                           col_points = rgb(0.5, 0.5, 0.5, 0.3),
                           col_gene_highlight = rgb(82, 185, 44, maxColorValue = 255),
                           cex_axis = 1.5, 
                           cex_lab = 1.5,
                           cex_points = 2.5,
                           density = 10,
                           gene_names = variable_selection_res$selected_variables,
                           lty_polygon = 1,
                           lwd_grid_major = 2,
                           lwd_grid_minor = 1,
                           lwd_axis = 1.5,
                           lwd_axis_ticks = 1.5,
                           lwd_polygon = 2,
                           lwd_polygon_bold = 4)
graphics.off()

# 
# round(variable_selection_res$logpval_vec[order(names(variable_selection_res$logpval_vec))],2)
# round(variable_selection_res$cor_vec_intial[order(names(variable_selection_res$cor_vec_intial))],2)
# idx <- order(variable_selection_res$logpval_vec, decreasing = T)[1:10]
# variable_selection_res$cor_vec_intial[idx]
