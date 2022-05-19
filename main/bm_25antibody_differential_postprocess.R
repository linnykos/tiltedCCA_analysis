rm(list=ls())
load("../../../out/main/citeseq_bm25_differential.RData")
load("../../../out/main/citeseq_bm25_tcca.RData")
source("bm_25antibody_colorPalette.R")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

# temporary fix
num_celltypes <- ceiling(sqrt(2*length(gene_de_list)))
combn_mat <- combn(num_celltypes, 2)
stopifnot(ncol(combn_mat) == length(gene_de_list))
gene_de_list <- list(de_list = gene_de_list, combn_mat = combn_mat, level_vec = 1:num_celltypes)
############

Seurat::DefaultAssay(bm) <- "RNA"
gene_names <- Seurat::VariableFeatures(bm)
logpval_vec <- sapply(1:length(gene_names), function(k){
  if(k %% floor(length(gene_names)/10) == 0) cat('*')
  gene <- gene_names[k]
  
  # cycle through all the celltypes
  celltype_vec <- sapply(1:length(gene_de_list$level_vec), function(i){
    idx <- which(gene_de_list$combn_mat == i, arr.ind = T)[,2]
    vec <-  sapply(idx, function(j){
      idx <- which(rownames(gene_de_list$de_list[[j]]) == gene)
      if(length(idx) == 0) return(1)
      gene_de_list$de_list[[j]][idx, "p_val"]
    })
    stats::quantile(vec, probs = 0.75)
  })
  
  max(-log10(celltype_vec))
})
names(logpval_vec) <- Seurat::VariableFeatures(bm)
logpval_vec <- pmin(logpval_vec, 300)

rsquare_vec <- tiltedCCA:::postprocess_alignment(input_obj = multiSVD_obj,
                                                 bool_use_denoised = T,
                                                 seurat_obj = bm,
                                                 input_assay = 1,
                                                 seurat_assay = "RNA",
                                                 seurat_slot = "data")
all(names(logpval_vec) == names(rsquare_vec))
stats::median(rsquare_vec[which(logpval_vec >= 1)])

png("../../../out/figures/main/citeseq_bm25_differential_gene.png",
    height = 3500, width = 2500, res = 500, units = "px")
# par(mar = c(5,5,4,1), bg = NA)
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "Human BM (CITE-Seq, RNA+ADT)\nGene differentiability vs. alignment",
                           bool_mark_ymedian = T,
                           col_gene_highlight_border = rgb(255, 205, 87, 255*0.5, maxColorValue = 255),
                           col_points = rgb(0.6, 0.6, 0.6, 0.1),
                           cex_axis = 1.5, 
                           cex_lab = 1.5,
                           cex_points = 2.5,
                           lty_polygon = 2,
                           lwd_grid_major = 2,
                           lwd_grid_minor = 1,
                           lwd_axis = 1.5,
                           lwd_axis_ticks = 1.5,
                           lwd_polygon_bold = 5,
                           mark_median_xthres = 1)
graphics.off()


png("../../../out/figures/main/citeseq_bm25_differential_gene_blank.png",
    height = 3500, width = 2500, res = 500, units = "px")
# par(mar = c(5,5,4,1), bg = NA)
par(mar = c(5,5,4,1))
tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                           logpval_vec = logpval_vec,
                           main = "Human BM (CITE-Seq, RNA+ADT)\nGene differentiability vs. alignment",
                           bool_hide_points = T,
                           bool_mark_ymedian = F,
                           cex_axis = 1.5, 
                           cex_lab = 1.5,
                           lwd_grid_major = 2,
                           lwd_grid_minor = 1,
                           lwd_axis = 1.5,
                           lwd_axis_ticks = 1.5)
graphics.off()

###################

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
bm$celltype.l3 <- celltype_l3
Seurat::Idents(bm) <- "celltype.l3"
level_vec <- sort(levels(Seurat::Idents(bm)))
num_celltype <- length(level_vec)
combn_mat <- utils::combn(num_celltype, 2)

p_val_thres <- 1e-4
log_thres <- 2
num_instances <- round(length(level_vec)*.25)
gene_list <- lapply(1:length(level_vec), function(i){
  idx <- which(combn_mat == i, arr.ind = T)[,2]
  lis <- lapply(gene_de_list[idx], function(mat){
    idx1 <- which(abs(mat[,2]) >= log_thres)
    idx2 <- which(abs(mat[,5]) <= p_val_thres)
    rownames(mat)[intersect(idx1, idx2)]
  })
  vec <- unlist(lis)
  tab <- table(vec)
  name_vec <- names(tab)[which(tab >= num_instances)]
  
  if(length(name_vec) <= 5){
    # find the genes that have the smallest median p-value
    all_genes <- unique(unlist(lapply(gene_de_list[idx], rownames)))
    all_genes <- all_genes[!all_genes %in% name_vec]
    
    p_val_vec <- sapply(all_genes, function(gene){
      val <- stats::median(sapply(idx, function(j){
        row_idx <- which(rownames(gene_de_list[[j]]) == gene)
        if(length(row_idx) == 0) return(1)
        gene_de_list[[j]][row_idx, "p_val"]
      }))
    })
    
    name_vec <- c(name_vec, all_genes[order(p_val_vec, decreasing = F)][1:(5-length(name_vec))])
  }
  
  name_vec
})
names(gene_list) <- level_vec
sapply(gene_list, length)

for(i in 1:length(gene_list)){
  png(paste0("../../../out/figures/main/citeseq_bm25_differential_gene_", names(gene_list)[i], ".png"),
      height = 3500, width = 2500, res = 500, units = "px")
  par(mar = c(5,5,4,1))
  tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                             logpval_vec = logpval_vec,
                             main = paste0("Human BM (CITE-Seq, RNA+ADT)\n", names(gene_list)[i], " genes"),
                             bool_mark_ymedian = F,
                             bool_polygon_mean = T,
                             col_points = rgb(0.5, 0.5, 0.5, 0.1),
                             col_gene_highlight = col_palette2[names(gene_list)[i]],
                             cex_axis = 1.5, 
                             cex_lab = 1.5,
                             cex_points = 2.5,
                             density = 10,
                             gene_names = gene_list[[i]],
                             lty_polygon = 1,
                             lwd_grid_major = 2,
                             lwd_grid_minor = 1,
                             lwd_axis = 1.5,
                             lwd_axis_ticks = 1.5,
                             lwd_polygon = 2,
                             lwd_polygon_bold = 4)
  graphics.off()
}

#######################

hk_genes <- read.csv("/home/stat/kevinl1/project/tiltedCCA/data/housekeeping.txt")[,1]
negative_list <- list(Housekeeping = hk_genes[which(hk_genes %in% names(logpval_vec))], 
                      Cell_cycle = c(cc.genes$s.genes[which(cc.genes$s.genes %in% names(logpval_vec))],
                                     cc.genes$g2m.genes[which(cc.genes$g2m.genes %in% names(logpval_vec))])
)

for(i in 1:length(negative_list)){
  png(paste0("../../../out/figures/main/citeseq_bm25_differential_gene_", names(negative_list)[i], ".png"),
      height = 3500, width = 2500, res = 500, units = "px")
  par(mar = c(5,5,4,1))
  tiltedCCA:::plot_alignment(rsquare_vec = rsquare_vec,
                             logpval_vec = logpval_vec,
                             main = paste0("Human BM (CITE-Seq, RNA+ADT)\n", names(negative_list)[i], " genes"),
                             bool_mark_ymedian = F,
                             bool_polygon_mean = T,
                             col_points = rgb(0.5, 0.5, 0.5, 0.1),
                             col_gene_highlight = "black",
                             cex_axis = 1.5, 
                             cex_lab = 1.5,
                             cex_points = 2.5,
                             density = 10,
                             gene_names = negative_list[[i]],
                             lty_polygon = 1,
                             lwd_grid_major = 2,
                             lwd_grid_minor = 1,
                             lwd_axis = 1.5,
                             lwd_axis_ticks = 1.5,
                             lwd_polygon = 2,
                             lwd_polygon_bold = 4)
  graphics.off()
}


