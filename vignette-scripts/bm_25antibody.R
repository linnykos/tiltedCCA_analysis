rm(list=ls())

# from what used to be https://satijalab.org/seurat/v4.0/weighted_nearest_neighbor_analysis.html
# now at https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-cite-seq-rna-adt-1

library(Seurat)
library(SeuratData)
library(cowplot)
library(dplyr)

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

# SeuratData::InstallData("bmcite")
bm <- SeuratData::LoadData(ds = "bmcite")

Seurat::DefaultAssay(bm) <- 'RNA'
bm <- Seurat::NormalizeData(bm)
bm <- Seurat::FindVariableFeatures(bm) 
bm <- Seurat::ScaleData(bm) 
bm <- Seurat::RunPCA(bm, 
                     verbose = F)

Seurat::DefaultAssay(bm) <- 'ADT'
Seurat::VariableFeatures(bm) <- rownames(bm[["ADT"]])
bm <- Seurat::NormalizeData(bm, 
                            normalization.method = 'CLR', 
                            margin = 2)
bm <- Seurat::ScaleData(bm) 
bm <- Seurat::RunPCA(bm, 
                     npcs = 18, 
                     reduction.name = 'apca', 
                     verbose = F)

set.seed(10)
bm <- Seurat::FindMultiModalNeighbors(
  bm, 
  reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:18), 
  modality.weight.name = "RNA.weight"
)

set.seed(10)
bm <- Seurat::RunUMAP(bm, 
                      reduction = 'pca', 
                      dims = 1:30, assay = 'RNA',
                      reduction.name = 'rna.umap', 
                      reduction.key = 'rnaUMAP_')
set.seed(10)
bm <- Seurat::RunUMAP(bm, 
                      reduction = 'apca', 
                      dims = 1:18, assay = 'ADT',
                      reduction.name = 'adt.umap', 
                      reduction.key = 'adtUMAP_')
set.seed(10)
bm <- Seurat::RunUMAP(bm,
                      nn.name = "weighted.nn", 
                      reduction.name = "wnn.umap", 
                      reduction.key = "wnnUMAP_")

################################

col_palette <- c(
  "CD14 Mono" = rgb(255, 108, 145, maxColorValue = 255),
  "CD16 Mono" = rgb(255, 104, 159, maxColorValue = 255),
  "CD4 Memory" = rgb(231, 134, 26, maxColorValue = 255),
  "CD4 Naive"= rgb(224, 139, 0, maxColorValue = 255),
  "CD56 bright NK" = rgb(0.7, 0.7, 0.7),
  "CD8 Effector_1" = rgb(2, 190, 108, maxColorValue = 255),
  "CD8 Effector_2" = rgb(0, 188, 89, maxColorValue = 255),
  "CD8 Memory_1" = rgb(57, 182, 0, maxColorValue = 255),
  "CD8 Memory_2" = rgb(0, 186, 66, maxColorValue = 255),
  "CD8 Naive" = rgb(1, 184, 31, maxColorValue = 255),
  "cDC2" = rgb(255, 98, 188, maxColorValue = 255),
  "gdT" = rgb(133, 173, 0, maxColorValue = 255),
  "GMP" = rgb(255, 97, 201, maxColorValue = 255),
  "HSC" = rgb(247, 99, 224, maxColorValue = 255),
  "LMPP" = rgb(240, 102, 234, maxColorValue = 255),
  "MAIT" = rgb(149, 169, 0, maxColorValue = 255),
  "Memory B" = rgb(0, 180, 239, maxColorValue = 255),
  "Naive B" = rgb(2, 165, 255, maxColorValue = 255),
  "NK" = rgb(0.6, 0.6, 0.6),
  "pDC" = rgb(255, 101, 174, maxColorValue = 255),
  "Plasmablast" = rgb(1, 191, 196, maxColorValue = 255),
  "Prog_B 1" = rgb(172, 136, 255, maxColorValue = 255),
  "Prog_B 2" = rgb(121, 151, 255, maxColorValue = 255),
  "Prog_DC" = rgb(252, 97, 213, maxColorValue = 255),
  "Prog_Mk" = rgb(231, 107, 243, maxColorValue = 255),
  "Prog_RBC" = rgb(220, 113, 250, maxColorValue = 255),
  "Treg" = rgb(243, 123, 89, maxColorValue = 255)
)

p1 <- Seurat::DimPlot(bm, 
                      reduction = "rna.umap",
                      group.by = "celltype.l2", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p1 <- p1 + Seurat::NoLegend()
p1 <- p1 + ggplot2::ggtitle("RNA") + ggplot2::labs(x = "", y = "")
p1 <- p1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p2 <- Seurat::DimPlot(bm, 
                      reduction = "adt.umap",
                      group.by = "celltype.l2", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p2 <- p2 + Seurat::NoLegend()
p2 <- p2 + ggplot2::ggtitle("ADT") + ggplot2::labs(x = "", y = "")
p2 <- p2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p3 <- Seurat::DimPlot(bm, 
                      reduction = "wnn.umap",
                      group.by = "celltype.l2", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p3 <- p3 + Seurat::NoLegend()
p3 <- p3 + ggplot2::ggtitle("WNN") + ggplot2::labs(x = "", y = "")
p3 <- p3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p_all <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
ggplot2::ggsave(filename = "vignettes/bm_wnn.png",
                p_all, device = "png", height = 1250, width = 3500, units = "px",
                dpi = 300)

############

Seurat::DefaultAssay(bm) <- "RNA"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:30)
bm <- Seurat::FindClusters(bm, resolution = 0.25)

Seurat::DefaultAssay(bm) <- "ADT"
set.seed(10)
bm <- Seurat::FindNeighbors(bm, dims = 1:18, reduction = "apca")
bm <- Seurat::FindClusters(bm, resolution = 0.25)

p1 <-Seurat::DimPlot(bm, 
                     reduction = "rna.umap",
                     group.by = "RNA_snn_res.0.25", 
                     label = TRUE, label.size = 2.5)
p1 <- p1 + Seurat::NoLegend() + ggplot2::labs(x = "", y = "")
p1 <- p1 + ggplot2::ggtitle("RNA clustering")
p1 <- p1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p2 <-Seurat::DimPlot(bm, 
                     reduction = "adt.umap",
                     group.by = "ADT_snn_res.0.25", 
                     label = TRUE, label.size = 2.5)
p2 <- p2 + Seurat::NoLegend() + ggplot2::labs(x = "", y = "")
p2 <- p2 + ggplot2::ggtitle("ADT clustering")
p2 <- p2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p_all <- cowplot::plot_grid(p1, p2, ncol = 2)
ggplot2::ggsave(filename = "vignettes/bm_clustering.png",
                p_all, device = "png", height = 1250, width = 2250, units = "px",
                dpi = 300)

#########

Seurat::DefaultAssay(bm) <- "RNA"
mat_1 <- Matrix::t(bm[["RNA"]]@data[Seurat::VariableFeatures(object = bm),])
Seurat::DefaultAssay(bm) <- "ADT"
mat_2 <- Matrix::t(bm[["ADT"]]@data)

mat_1b <- mat_1
sd_vec <- sparseMatrixStats::colSds(mat_1b)
if(any(sd_vec <= 1e-6)){
  mat_1b <- mat_1b[,-which(sd_vec <= 1e-6)]
}

mat_2b <- mat_2
sd_vec <- sparseMatrixStats::colSds(mat_2b)
if(any(sd_vec <= 1e-6)){
  mat_2b <- mat_2b[,-which(sd_vec <= 1e-6)]
}

set.seed(10)
multiSVD_obj <- tiltedCCA:::create_multiSVD(mat_1 = mat_1b, mat_2 = mat_2b,
                                            dims_1 = 1:30, dims_2 = 1:18,
                                            center_1 = T, center_2 = T,
                                            normalize_row = T,
                                            normalize_singular_value = T,
                                            recenter_1 = F, recenter_2 = F,
                                            rescale_1 = F, rescale_2 = F,
                                            scale_1 = T, scale_2 = T,
                                            verbose = 1)
multiSVD_obj <- tiltedCCA:::form_metacells(input_obj = multiSVD_obj,
                                           large_clustering_1 = as.factor(bm$RNA_snn_res.0.25), 
                                           large_clustering_2 = as.factor(bm$ADT_snn_res.0.25), 
                                           num_metacells = 5000,
                                           verbose = 1)
multiSVD_obj <- tiltedCCA:::compute_snns(input_obj = multiSVD_obj,
                                         latent_k = 20,
                                         num_neigh = 30,
                                         bool_cosine = T,
                                         bool_intersect = F,
                                         min_deg = 30,
                                         verbose = 2)
multiSVD_obj <- tiltedCCA:::tiltedCCA(input_obj = multiSVD_obj,
                                      verbose = 1)
multiSVD_obj <- tiltedCCA:::fine_tuning(input_obj = multiSVD_obj,
                                        verbose = 1)
# takes 20 minutes
multiSVD_obj <- tiltedCCA:::tiltedCCA_decomposition(input_obj = multiSVD_obj,
                                                    verbose = 1)


set.seed(10)
bm[["common_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                    what = "common",
                                                    aligned_umap_assay = "rna.umap",
                                                    seurat_obj = bm,
                                                    seurat_assay = "RNA",
                                                    verbose = 1)
set.seed(10)
bm[["distinct1_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                       what = "distinct_1",
                                                       aligned_umap_assay = "rna.umap",
                                                       seurat_obj = bm,
                                                       seurat_assay = "RNA",
                                                       verbose = 1)
set.seed(10)
bm[["distinct2_tcca"]] <- tiltedCCA:::create_SeuratDim(input_obj = multiSVD_obj,
                                                       what = "distinct_2",
                                                       aligned_umap_assay = "rna.umap",
                                                       seurat_obj = bm,
                                                       seurat_assay = "RNA",
                                                       verbose = 1)

p1 <- Seurat::DimPlot(bm, 
                      reduction = "common_tcca",
                      group.by = "celltype.l2", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p1 <- p1 + Seurat::NoLegend()
p1 <- p1 + ggplot2::ggtitle("Common") + ggplot2::labs(x = "", y = "")
p1 <- p1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p2 <- Seurat::DimPlot(bm, 
                      reduction = "distinct1_tcca",
                      group.by = "celltype.l2", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p2 <- p2 + Seurat::NoLegend()
p2 <- p2 + ggplot2::ggtitle("RNA Distinct") + ggplot2::labs(x = "", y = "")
p2 <- p2 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p3 <- Seurat::DimPlot(bm, 
                      reduction = "distinct2_tcca",
                      group.by = "celltype.l2", 
                      cols = col_palette, 
                      label = T, repel = T, label.size = 2.5)
p3 <- p3 + Seurat::NoLegend()
p3 <- p3 + ggplot2::ggtitle("ADT Distinct") + ggplot2::labs(x = "", y = "")
p3 <- p3 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))

p_all <- cowplot::plot_grid(p1, p2, p3, ncol = 3)
ggplot2::ggsave(filename = "vignettes/bm_tiltedcca.png",
                p_all, device = "png", height = 1250, width = 3500, units = "px",
                dpi = 300)

#######################

celltype_l3 <- as.character(bm$celltype.l2)
celltype_l3[celltype_l3 %in% c("NK", "CD56 bright NK")] <- "NK_all"
celltype_l3[celltype_l3 %in% c("MAIT", "gdT")] <- "MAIT-gdT"
celltype_l3[celltype_l3 %in% c("Memory B", "Naive B")] <- "B"
celltype_l3[celltype_l3 %in% c("Prog_B 1", "Prog_B 2")] <- "Prog_B"
celltype_l3[celltype_l3 %in% c("Prog_Mk", "Plasmablast", "LMPP", "Treg")] <- NA
col_palette2 <- c(col_palette, 
                  col_palette["NK"], 
                  col_palette["MAIT"], 
                  col_palette["Memory B"], 
                  col_palette["Prog_B 1"])
names(col_palette2) <- c(names(col_palette), "NK_all", "MAIT-gdT", "B", "Prog_B")

keep_celltypes <- c("CD4 Naive", "CD8 Naive", "B", "CD14 Mono")
celltype_l3[!celltype_l3 %in% keep_celltypes] <- NA

bm$celltype.l3 <- celltype_l3
# start: 10:59
print("Working on RNA")
gene_de_list <- tiltedCCA:::differential_expression(seurat_obj = bm,
                                                    assay = "RNA",
                                                    idents = "celltype.l3",
                                                    test_use = "MAST",
                                                    slot = "counts")

logpval_vec <- postprocess_depvalue(de_list = gene_de_list, 
                                    maximum_output = 300)
rsquare_vec <- postprocess_modality_alignment(input_obj = multiSVD_obj,
                                              bool_use_denoised = T,
                                              seurat_obj = bm,
                                              input_assay = 1,
                                              seurat_assay = "RNA",
                                              seurat_slot = "data")

col_palette_enrichment <- paste0(grDevices::colorRampPalette(c('lightgrey', 'blue'))(100), "33")
gene_breaks <- seq(0, 0.18, length.out = 100)
gene_depth <- Matrix::rowSums(bm[["RNA"]]@counts[names(rsquare_vec),])/ncol(bm)
gene_depth <- log10(gene_depth+1)
gene_color <- sapply(gene_depth, function(x){
  col_palette_enrichment[which.min(abs(gene_breaks - x))]
})
ord_idx <- order(gene_depth, decreasing = F)

gene_list <- postprocess_marker_variables(de_list)
cycling_genes <- c(cc.genes$s.genes[which(cc.genes$s.genes %in% names(logpval_vec))],
                   cc.genes$g2m.genes[which(cc.genes$g2m.genes %in% names(logpval_vec))])

png("vignettes/bm_gene-alignment.png", height = 1250, width = 3500,
    units = "px", res = 300)
par(mfrow = c(1,4), mar = c(5,5,5,1))
plot_alignment(rsquare_vec = rsquare_vec[ord_idx],
               logpval_vec = logpval_vec[ord_idx],
               main = "Human BM (CITE-Seq, RNA+ADT)\nGene differentiability vs. alignment",
               bool_mark_ymedian = T,
               col_gene_highlight_border = rgb(255, 205, 87, 255*0.5, maxColorValue = 255),
               col_points = gene_color[ord_idx],
               mark_median_xthres = 10,
               lwd_polygon_bold = 3)

for(i in 1:2){
  plot_alignment(rsquare_vec = rsquare_vec,
                 logpval_vec = logpval_vec,
                 main = paste0("Human BM (CITE-Seq, RNA+ADT)\n", names(gene_list)[i], " genes"),
                 col_gene_highlight = col_palette2[names(gene_list)[i]],
                 gene_names = gene_list[[i]],
                 lwd_polygon_bold = 3)
}

plot_alignment(rsquare_vec = rsquare_vec,
               logpval_vec = logpval_vec,
               main = paste0("Human BM (CITE-Seq, RNA+ADT)\nCell cycling genes"),
               col_gene_highlight = "black",
               gene_names = cycling_genes,
               lwd_polygon_bold = 3)
graphics.off()

###################

membership_vec <- factor(bm$celltype.l2)

# start 11:27. takes 5 minutes
cell_enrichment_res <- tiltedCCA::postprocess_cell_enrichment(
  input_obj = multiSVD_obj, 
  membership_vec = membership_vec, 
  max_subsample = 1000,
  verbose = 1
)

png("vignettes/bm_cell_enrichment_common-vs-distinct.png",
    height = 2500, width = 5000, res = 500, units = "px")
plot_cell_enrichment(cell_enrichment_res,
                     cex_axis = 1.25,
                     cex_lab = 1.5,
                     xlab_1 = "RNA Distinct",
                     xlab_2 = "ADT Distinct")
graphics.off()

#############

# saving
multiSVD_obj$common_mat_1 <- NULL
multiSVD_obj$common_mat_2 <- NULL
multiSVD_obj$distinct_mat_1 <- NULL
multiSVD_obj$distinct_mat_2 <- NULL

save(multiSVD_obj, file = "../../out/experiment/vignette/bm_25antibody.RData")

save(gene_de_list, adt_de_list,
     file = "../../out/experiment/vignette/bm_25antibody_de.RData")

