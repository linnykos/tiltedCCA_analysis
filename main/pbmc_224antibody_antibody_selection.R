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
  logpval_vec = logpval_vec,
  cor_threshold = 0.45,
  input_assay = 2,
  max_variables = 10,
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
pbmc[["ADT2"]]@var.features <- variable_selection_res$selected_variables

Seurat::DefaultAssay(pbmc) <- "ADT2"
pbmc <- Seurat::RunPCA(pbmc, reduction.name = 'apca2', 
                       npcs = length(variable_selection_res$selected_variables), 
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

save(variable_selection_res, pbmc, logpval_vec,
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

################################

multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(multiSVD_obj, 
                                                    bool_modality_1_full = F,
                                                    bool_modality_2_full = F,
                                                    verbose = 0)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
reference_dimred <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "common_dimred")
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
distinct_mat <- tiltedCCA:::.get_tCCAobj(multiSVD_obj, apply_postDimred = F, what = "distinct_mat")

rsquare_vec <- sapply(1:ncol(distinct_mat), function(j){
  tiltedCCA:::.linear_regression(bool_include_intercept = T,
                                 bool_center_x = T,
                                 bool_center_y = T,
                                 bool_scale_x = T,
                                 bool_scale_y = T,
                                 return_type = "r_squared", 
                                 x_mat = reference_dimred,
                                 y_vec = distinct_mat[,j])
})
names(rsquare_vec) <- colnames(distinct_mat)
rsquare_vec[variable_selection_res$selected_variables]

# rsquare_vec <- tiltedCCA:::postprocess_alignment(input_obj = multiSVD_obj,
#                                                  bool_use_denoised = T,
#                                                  seurat_obj = pbmc,
#                                                  input_assay = 2,
#                                                  min_subsample_cell = NULL,
#                                                  seurat_assay = "SCT",
#                                                  seurat_celltype_variable = "celltype.l2",
#                                                  seurat_slot = "data",
#                                                  verbose = 2)
# all(names(logpval_vec) == names(rsquare_vec))

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


