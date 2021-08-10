rm(list=ls()); set.seed(10)

library(Seurat); library(multiomicCCA)
# library(Signac); library(EnsDb.Hsapiens.v86); 
# library(BSgenome.Hsapiens.UCSC.hg38); library(GenomeInfoDb)
# library(dplyr); library(ggplot2)

load("../../../../out/Writeup10/Writeup10_10x_pbmc_preprocess4.RData")
date_of_run <- Sys.time(); session_info <- sessionInfo()

load("../../../../out/Writeup14b/Writeup14b_10x_pbmc_dcca_laplacian_variablecalculations.RData")

# marker_genes <- c("AVP", "LMO4", "PF4", "BLVRB", "MME", "DERL3", "CLEC9A", 
#                   "CDC1", "MPO", "AZU1", "CD14", "FCGR3A", "VREB3", "MS4A1", 
#                   "CD79A", "IGKC", "PF4", "XCL1", "CD8A", "CD4", "SH2D1A")
# marker_idx <- which(marker_genes %in% rownames(pbmc[["RNA"]]))
# # see if we can fix the unmatched ones
# marker_unmatched <- setdiff(1:length(marker_genes), marker_idx)
# gene_syn <- Seurat::GeneSymbolThesarus(marker_genes[marker_unmatched], several.ok = T)
# marker_idx2 <- lapply(gene_syn, function(x){
#   vec <- which(x %in% rownames(pbmc[["RNA"]]))
#   if(length(vec) > 0){
#     return(x[vec[1]])
#   } else {
#     return(numeric(0))
#   }
# })
# marker_genes <- c(marker_genes[marker_idx], unlist(marker_idx2))
# names(marker_genes) <- NULL
# 
# # now find which of the genes have been selected by the 3000 genes 
# marker_genes <- marker_genes[which(marker_genes %in% rownames(dcca_res$svd_1$v))]

#####################


val_vec <- sapply(gene_smoothed, function(x){
  (x$d_variance - x$c_variance)/x$e_variance
})
name_vec <- rownames(dcca_res$svd_1$v)
p1 <- length(name_vec)
factor_vec <- rep(0, p1); 
factor_vec[which(sapply(gene_smoothed, function(x){min(x$c_r2, x$d_r2) > 0.6}))] <- 2
idx <- which(factor_vec == 2)
idx2 <- which.max(val_vec[idx])
factor_vec[idx[idx2]] <- 1; factor_vec <- as.factor(factor_vec)
col_vec <- c("black", "red", "green"); names(col_vec) <- c("0", "1", "2")
tmp <- multiomicCCA::plot_laplacian_variables(val_vec, name_vec, factor_vec, col_vec,
                                              ylab = "Distinct-common (normalized)", main = "RNA enrichment",
                                              text_cex = 4)
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_rna_var_enrichment.png",
                tmp, device = "png", width = 3, height = 3, units = "in")

################

val_vec <- sapply(atac_smoothed, function(x){
  (x$d_variance - x$c_variance)/x$e_variance
})
name_vec <- rownames(dcca_res$svd_2$v)
p2 <- length(name_vec)
factor_vec <- rep(0, p2); 
factor_vec[which(sapply(atac_smoothed, function(x){min(x$c_r2, x$d_r2) > 0.6}))] <- 2
idx <- which(factor_vec == 2)
idx2 <- which.max(val_vec[idx])
factor_vec[idx[idx2]] <- 1; factor_vec <- as.factor(factor_vec)
col_vec <- c("black", "red", "green"); names(col_vec) <- c("0", "1", "2")
tmp <- multiomicCCA::plot_laplacian_variables(val_vec, name_vec, factor_vec, col_vec,
                                              ylab = "Distinct-common (normalized)", main = "ATAC enrichment",
                                              text_cex = 4)
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14b/Writeup14b_10x_pbmc_atac_var_enrichment.png",
                tmp, device = "png", width = 3, height = 3, units = "in")







