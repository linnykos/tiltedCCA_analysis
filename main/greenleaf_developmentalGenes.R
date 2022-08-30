rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_greenleaf_tcca_RNA-geneActivity.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

colnames(multiSVD_obj$common_mat_2) <- sapply(colnames(multiSVD_obj$common_mat_2), function(x){
  strsplit(x, split = "ATAC-")[[1]][2]
})
colnames(multiSVD_obj$distinct_mat_2) <- sapply(colnames(multiSVD_obj$distinct_mat_2), function(x){
  strsplit(x, split = "ATAC-")[[1]][2]
})

cell_idx <- unique(c(which(greenleaf$Lineage1 == 1), which(greenleaf$Lineage2 == 1)))

potential_genes <- sort(intersect(colnames(multiSVD_obj$common_mat_1), 
                                  colnames(multiSVD_obj$common_mat_2)))
bool_vec <- sapply(1:length(potential_genes), function(i){
  if(i %% floor(length(potential_genes)/10) == 0) cat('*')
  gene_name <- potential_genes[i]
  quantile(as.numeric(greenleaf[["SCT"]]@data[gene_name,cell_idx]), probs = 0.9) > 0
})
variable_names <- potential_genes[which(bool_vec)]

selection_res <- tiltedCCA:::postprocess_smooth_variable_selection(
  input_obj = multiSVD_obj,
  bool_use_denoised = T,
  bool_include_intercept = T,
  bool_use_metacells = F,
  bool_use_both_modalities = T,
  cell_idx = cell_idx,
  cor_threshold = 0.95,
  num_variables = 50,
  sd_quantile = 0.05,
  seurat_obj = greenleaf,
  seurat_assay_1 = "SCT",
  seurat_assay_2 = "customGAct",
  seurat_slot = "data",
  variable_names = variable_names,
  verbose = 2
)
sort(selection_res$selected_variables)
sort(names(selection_res$alignment_1))

Seurat::DefaultAssay(greenleaf) <- "SCT"
plot1 <- Seurat::FeaturePlot(greenleaf, 
                             features = sort(selection_res$selected_variables)[1:25],
                             reduction = "common_tcca",
                             ncol = 5)
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_develompentalGenes1.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")

plot1 <- Seurat::FeaturePlot(greenleaf, 
                             features = sort(selection_res$selected_variables)[26:50],
                             reduction = "common_tcca",
                             ncol = 5)
ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_tcca_develompentalGenes2.png"),
                plot1, device = "png", width = 20, height = 15, units = "in")

save(selection_res, date_of_run, session_info, cell_idx,
     file = "../../../out/main/10x_greenleaf_developmentalGenes.RData")

sink("../../../out/main/10x_greenleaf_developmentalGenes.txt")
for(i in 1:length(selection_res$selected_variables)){
  cat(selection_res$selected_variables[i])
  cat("\n")
}
sink()