rm(list=ls())
load("../../../../out/Writeup14b/Writeup14b_mouseicb_dcca_tmp.RData")
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
idx4 <- which(myeloid2@meta.data$celltype == "R499_dICB")
start_idx <- intersect(intersect(idx1, intersect(idx2, idx3)), idx4)
table(myeloid2@meta.data$celltype[start_idx])

multiomicCCA::compute_enrichment_scores(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx)

set.seed(10)
res <- multiomicCCA::detect_rare_cell(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx,
                                      common_enrich = F, distinct_enrich_1 = T, 
                                      distinct_enrich_2 = F, 
                                      custom_threshold = c(0.25, 0.4, 0.1),
                                      deg_threshold = 0, max_size = 400,
                                      verbose = T)

multiomicCCA::compute_enrichment_scores(combined_g, rna_frnn$d_g, atac_frnn$d_g, res$idx)
table(myeloid2@meta.data$celltype[res$idx])

####################

full_col_vec <- scales::hue_pal()(4)[as.numeric(as.factor(myeloid2@meta.data$celltype))]
n <- length(full_col_vec)

png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct1.png", 
    width = 3000, height = 2000, units = "px", res = 300)
.custom_plot1(myeloid2, full_col_vec, res)
graphics.off()


png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct1b.png", 
    width = 3000, height = 2000, units = "px", res = 300)
.custom_plot2(myeloid2, full_col_vec, res)
graphics.off()

##############################
##############################
##############################

tmp <- myeloid2[["distinct"]]@cell.embeddings
# idx1 <- which(myeloid2@meta.data$celltype == "B16_dICB")
idx1 <- which(tmp[,1] <= -3.5)
idx2 <- which(tmp[,2] <= -3.5)
idx3 <- which(myeloid2@meta.data$celltype == "R499_dICB")
start_idx <- intersect(intersect(idx1, idx2), idx3)
table(myeloid2@meta.data$celltype[start_idx])

multiomicCCA::compute_enrichment_scores(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx)

set.seed(10)
res <- multiomicCCA::detect_rare_cell(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx,
                                      common_enrich = F, distinct_enrich_1 = T, 
                                      distinct_enrich_2 = F, 
                                      deg_threshold = 0, max_size = 500,
                                      verbose = T)
table(myeloid2@meta.data$celltype[res$idx])
length(res$idx)/length(res$baseline_idx)

full_col_vec <- scales::hue_pal()(4)[as.numeric(as.factor(myeloid2@meta.data$celltype))]
n <- length(full_col_vec)

png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct2.png", 
    width = 3000, height = 2000, units = "px", res = 300)
.custom_plot1(myeloid2, full_col_vec, res)
graphics.off()


png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct2b.png", 
    width = 3000, height = 2000, units = "px", res = 300)
.custom_plot2(myeloid2, full_col_vec, res)
graphics.off()


##############################
##############################
##############################

tmp <- myeloid2[["distinct2"]]@cell.embeddings
idx1 <- which(tmp[,2] >= 3.2)
idx2 <- which(myeloid2@meta.data$celltype == "B16_dICB")
start_idx <- intersect(idx1, idx2)
table(myeloid2@meta.data$celltype[start_idx])

multiomicCCA::compute_enrichment_scores(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx)

set.seed(10)
res <- multiomicCCA::detect_rare_cell(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx,
                                      common_enrich = F, distinct_enrich_1 = F, 
                                      distinct_enrich_2 = T, 
                                      deg_threshold = 0, max_size = 250,
                                      verbose = T)
table(myeloid2@meta.data$celltype[res$idx])
length(res$idx)/length(res$baseline_idx)

full_col_vec <- scales::hue_pal()(4)[as.numeric(as.factor(myeloid2@meta.data$celltype))]
n <- length(full_col_vec)

png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct3.png", 
    width = 3000, height = 2000, units = "px", res = 300)
.custom_plot1(myeloid2, full_col_vec, res)
graphics.off()


png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct3b.png", 
    width = 3000, height = 2000, units = "px", res = 300)
.custom_plot2(myeloid2, full_col_vec, res)
graphics.off()

#######################


tmp <- myeloid2[["combined"]]@cell.embeddings
idx1 <- which(tmp[,1] <= -3)
idx2 <- which(tmp[,1] >= -3.5)
idx3 <- which(tmp[,2] >= 0)
idx4 <- which(tmp[,2] <= 0.5)
idx5 <- which(myeloid2@meta.data$celltype == "R499_dICB")
start_idx <- Reduce(intersect, list(idx1, idx2, idx3, idx4, idx5))
table(myeloid2@meta.data$celltype[start_idx])

multiomicCCA::compute_enrichment_scores(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx)

set.seed(10)
res <- multiomicCCA::detect_rare_cell(combined_g, rna_frnn$d_g, atac_frnn$d_g, start_idx,
                                      common_enrich = T, distinct_enrich_1 = F, 
                                      distinct_enrich_2 = F, 
                                      custom_threshold = c(0.2, 0.1, 0.05),
                                      deg_threshold = 0, max_size = 250,
                                      verbose = T)
table(myeloid2@meta.data$celltype[res$idx])
length(res$idx)/length(res$baseline_idx)


full_col_vec <- scales::hue_pal()(4)[as.numeric(as.factor(myeloid2@meta.data$celltype))]
n <- length(full_col_vec)

png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_common1.png", 
    width = 3000, height = 2000, units = "px", res = 300)
.custom_plot1(myeloid2, full_col_vec, res)
graphics.off()


png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_common1b.png", 
    width = 3000, height = 2000, units = "px", res = 300)
.custom_plot2(myeloid2, full_col_vec, res)
graphics.off()

