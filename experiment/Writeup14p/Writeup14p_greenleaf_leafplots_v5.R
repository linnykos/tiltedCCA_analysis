rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/Writeup14p/Writeup14p_10x_greenleaf_tcca.RData")

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_common <- multiSVD_obj$common_mat_1[cell_idx[order_vec],]

load("../../../../out/main/10x_greenleaf_developmentalGenes.RData")
mentioned_genes <- c("ALDH2", "APOE", "AQP4", "ASCL1", "BHLHE22",
                     "BHLHE40", "C16orf89", "CAV2", "DLX2", "DOK5",
                     "DUSP1", "EOMES", "ETV4", "FOS", "FOXJ1", "GLI3",
                     "HAS2", "HES1", "HES4", "HSPA1A", "HSPA1B",
                     "ID3", "IGFBP7", "JUN", "KIF1A", "LIMCH1",
                     "MBP", "MEF2C", "NEUROD1", "NEUROD2", "NEUROD4",
                     "NEUROD6", "NEUROG1", "NEUROG2", "NFIA", "NFIB", "NFIC",
                     "NHLH1", "NR2F1", "PAX6", "RFX4", "RUNX1",
                     "OLIG1", "OLIG2", "SOX2", "SOX3", "SOX6",
                     "SOX9", "SOX10", "SOX21", "SPARCL1", "SNCB", "TBX",
                     "TNC", "TOP2A", "TRB1", "WNT11",
                     "SATB2")
mentioned_genes <- unique(c(mentioned_genes, selection_res$selected_variables))
mentioned_genes <- intersect(mentioned_genes, colnames(multiSVD_obj$common_mat_1))
mentioned_genes <- sort(unique(mentioned_genes))
rna_common <- rna_common[,mentioned_genes]

############

rna_common <- multiSVD_obj$common_mat_1
rna_distinct <- multiSVD_obj$distinct_mat_1
rna_mat <- multiSVD_obj$common_mat_1 + multiSVD_obj$distinct_mat_1
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 1)
svd_1 <- tiltedCCA:::.get_SVD(multiSVD_obj)
multiSVD_obj <- tiltedCCA:::.set_defaultAssay(multiSVD_obj, assay = 2)
svd_2 <- tiltedCCA:::.get_SVD(multiSVD_obj)
tmp <- crossprod(svd_2$u, svd_1$u)
svd_tmp <- svd(tmp)
rotation_mat <- tcrossprod(svd_tmp$u, svd_tmp$v)
atac_mat_full <- tcrossprod(tiltedCCA:::.mult_mat_vec(svd_2$u %*% rotation_mat, svd_1$d), svd_1$v)
n <- nrow(atac_mat_full)
atac_mat_common <- multiSVD_obj$tcca_obj$common_score %*% crossprod(multiSVD_obj$cca_obj$score_2, atac_mat_full)/n
atac_mat_distinct <- multiSVD_obj$tcca_obj$distinct_score_2 %*% crossprod(multiSVD_obj$cca_obj$score_2, atac_mat_full)/n

cell_idx <- intersect(which(greenleaf$Lineage1 == 1), which(!is.na(greenleaf$pseudotime)))
pseudotime_vec <- greenleaf$pseudotime[cell_idx]
order_vec <- order(pseudotime_vec, decreasing = F)
rna_mat <- rna_mat[cell_idx[order_vec],mentioned_genes]
rna_common <- rna_common[cell_idx[order_vec],mentioned_genes]
rna_distinct <- rna_distinct[cell_idx[order_vec],mentioned_genes]
atac_mat_full <- atac_mat_full[cell_idx[order_vec],mentioned_genes]
atac_mat_common <- atac_mat_common[cell_idx[order_vec],mentioned_genes]
atac_mat_distinct <- atac_mat_distinct[cell_idx[order_vec],mentioned_genes]

save(rna_mat, rna_common, rna_distinct, 
     atac_mat_full, atac_mat_common, atac_mat_distinct,
     file = "../../../../out/Writeup14p/tmp2.RData")

#################

load("../../out/experiment/Writeup14p/tmp2.RData")
n <- nrow(rna_mat)

plot(rna_mat[,"EOMES"], pch = 16)
points(atac_mat_common[,"EOMES"], pch = 16, col = 2)

prop <- abs(atac_mat_common[,"EOMES"])/(abs(atac_mat_common[,"EOMES"]) + abs(atac_mat_distinct[,"EOMES"]))
atac_vec <- atac_mat_full[,"EOMES"]*prop
plot(rna_mat[,"EOMES"], pch = 16)
points(atac_vec, pch = 16, col = 2)

# prop <- abs(atac_mat_common[,"EOMES"])/abs(atac_mat_distinct[,"EOMES"])
# atac_vec <- atac_mat_full[,"EOMES"]*prop
atac_vec <- rna_distinct[,"EOMES"]
rna_vec <- rna_mat[,"EOMES"]
atac_vec <- pmin(pmax(atac_vec, -5),10)
plot(atac_vec, pch = 16, col = 2, ylim = quantile(c(atac_vec, rna_vec), probs = c(0.01,0.99)))
points(rna_vec, pch = 16, col = 1)

base_palette <- RColorBrewer::brewer.pal(11, name = "RdYlBu")
color_vec <- rev(grDevices::colorRampPalette(base_palette)(n))
plot(atac_vec, rna_vec, col = color_vec, pch = 16)




