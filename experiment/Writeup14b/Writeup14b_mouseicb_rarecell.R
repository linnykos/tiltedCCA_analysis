rm(list=ls())
load("../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")
source("custom_plotting.R")

# nn <- 30
# set.seed(10)
# combined_g <- multiomicCCA::combine_frnn(dcca_res, g_1 = rna_frnn$c_g,
#                                          g_2 = atac_frnn$c_g, nn = nn, 
#                                          verbose = 2)
# rm(list = "nn")
# save.image("../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

tmp <- myeloid2[["distinct"]]@cell.embeddings
# idx1 <- which(myeloid2@meta.data$celltype == "B16_dICB")
idx1 <- which(tmp[,1] <= -4)
idx2 <- which(tmp[,2] >= -2.5)
idx3 <- which(tmp[,2] <= -1.5)
start_idx <- intersect(idx1, intersect(idx2, idx3))
table(myeloid2@meta.data$celltype[start_idx])

multiomicCCA::compute_enrichment_scores(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx)

set.seed(10)
res <- multiomicCCA::detect_rare_cell(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx,
                                      common_enrich = F, distinct_enrich_1 = T, 
                                      distinct_enrich_2 = F, 
                                      deg_threshold = 0, max_size = NA,
                                      verbose = T)

multiomicCCA::compute_enrichment_scores(combined_g, rna_frnn$d_g, atac_frnn$d_g, res$idx)
table(myeloid2@meta.data$celltype[res$idx])

####################

full_col_vec <- scales::hue_pal()(4)[as.numeric(as.factor(myeloid2@meta.data$celltype))]
n <- length(full_col_vec)

png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct1.png", 
    width = 3000, height = 2000, units = "px", res = 300)
graphics.off()


png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct1b.png", 
    width = 3000, height = 2000, units = "px", res = 300)
graphics.off()

##############################
##############################
##############################

tmp <- myeloid2[["distinct"]]@cell.embeddings
# idx1 <- which(myeloid2@meta.data$celltype == "B16_dICB")
idx1 <- which(tmp[,1] <= -3.5)
idx2 <- which(tmp[,2] <= -3.5)
start_idx <- intersect(idx1, idx2)
table(myeloid2@meta.data$celltype[start_idx])

multiomicCCA::compute_enrichment_scores(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx)

set.seed(10)
res <- multiomicCCA::detect_rare_cell(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx,
                                      common_enrich = F, distinct_enrich_1 = T, 
                                      distinct_enrich_2 = F, 
                                      custom_threshold = c(0.15, 0.9, 0.05),
                                      deg_threshold = 0.05, max_size = 500,
                                      verbose = T)
table(myeloid2@meta.data$celltype[res$idx])
length(res$idx)/length(res$baseline_idx)

full_col_vec <- scales::hue_pal()(4)[as.numeric(as.factor(myeloid2@meta.data$celltype))]
n <- length(full_col_vec)

png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct2.png", 
    width = 3000, height = 2000, units = "px", res = 300)
par(mfcol = c(2,3), mar = c(4,4,4,0.1))
layout <- myeloid2[["combined"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, 
           main = paste0("Common combined,\nInitial (", round(res$baseline_scores[1], 2), ")"))
.plot_umap(layout, full_col_vec, res$idx, 
           main = paste0("Common combined,\nFinal (", round(res$scores[1], 2), ")"))

layout <- myeloid2[["distinct"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, 
           main = paste0("RNA distinct,\nInitial (", round(res$baseline_scores[2], 2), ")"))
.plot_umap(layout, full_col_vec, res$idx, 
           main = paste0("RNA distinct,\nFinal (", round(res$scores[2], 2), ")"))

layout <- myeloid2[["distinct2"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx,
           main = paste0("ATAC distinct,\nInitial (", round(res$baseline_scores[3], 2), ")"))
.plot_umap(layout, full_col_vec, res$idx,
           main = paste0("ATAC distinct,\nFinal (", round(res$scores[3], 2), ")"))
graphics.off()


png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct2b.png", 
    width = 3000, height = 2000, units = "px", res = 300)
par(mfcol = c(2,3), mar = c(4,4,4,0.1))
layout <- myeloid2[["everything"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, "RNA everything,\nInitial")
.plot_umap(layout, full_col_vec, res$idx, "RNA everything,\nFinal")

layout <- myeloid2[["everything2"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, "ATAC everything,\nInitial")
.plot_umap(layout, full_col_vec, res$idx, "ATAC everything,\nFinal")

layout <- myeloid2[["both"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, "Both everything,\nInitial")
.plot_umap(layout, full_col_vec, res$idx, "Both everything,\nFinal")
graphics.off()


##############################
##############################
##############################

tmp <- myeloid2[["distinct2"]]@cell.embeddings
# idx1 <- which(myeloid2@meta.data$celltype == "B16_dICB")
start_idx <- which(tmp[,2] >= 3.2)
table(myeloid2@meta.data$celltype[start_idx])

multiomicCCA::compute_enrichment_scores(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx)

set.seed(10)
res <- multiomicCCA::detect_rare_cell(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx,
                                      common_enrich = F, distinct_enrich_1 = F, 
                                      distinct_enrich_2 = T, 
                                      custom_threshold = c(0.1, 0.05, 0.5),
                                      deg_threshold = 0, max_size = 250,
                                      verbose = T)
table(myeloid2@meta.data$celltype[res$idx])
length(res$idx)/length(res$baseline_idx)

full_col_vec <- scales::hue_pal()(4)[as.numeric(as.factor(myeloid2@meta.data$celltype))]
n <- length(full_col_vec)

png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct3.png", 
    width = 3000, height = 2000, units = "px", res = 300)
par(mfcol = c(2,3), mar = c(4,4,4,0.1))
layout <- myeloid2[["combined"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, 
           main = paste0("Common combined,\nInitial (", round(res$baseline_scores[1], 2), ")"))
.plot_umap(layout, full_col_vec, res$idx, 
           main = paste0("Common combined,\nFinal (", round(res$scores[1], 2), ")"))

layout <- myeloid2[["distinct"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, 
           main = paste0("RNA distinct,\nInitial (", round(res$baseline_scores[2], 2), ")"))
.plot_umap(layout, full_col_vec, res$idx, 
           main = paste0("RNA distinct,\nFinal (", round(res$scores[2], 2), ")"))

layout <- myeloid2[["distinct2"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx,
           main = paste0("ATAC distinct,\nInitial (", round(res$baseline_scores[3], 2), ")"))
.plot_umap(layout, full_col_vec, res$idx,
           main = paste0("ATAC distinct,\nFinal (", round(res$scores[3], 2), ")"))
graphics.off()


png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct3b.png", 
    width = 3000, height = 2000, units = "px", res = 300)
par(mfcol = c(2,3), mar = c(4,4,4,0.1))
layout <- myeloid2[["everything"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, "RNA everything,\nInitial")
.plot_umap(layout, full_col_vec, res$idx, "RNA everything,\nFinal")

layout <- myeloid2[["everything2"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, "ATAC everything,\nInitial")
.plot_umap(layout, full_col_vec, res$idx, "ATAC everything,\nFinal")

layout <- myeloid2[["both"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, "Both everything,\nInitial")
.plot_umap(layout, full_col_vec, res$idx, "Both everything,\nFinal")
graphics.off()


