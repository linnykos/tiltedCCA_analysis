rm(list=ls())
load("../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")

nn <- 30
set.seed(10)
combined_g <- multiomicCCA::combine_frnn(dcca_res, g_1 = rna_frnn$c_g,
                                         g_2 = atac_frnn$c_g, nn = nn, 
                                         verbose = 2)

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
                                      distinct_enrich_2 = F, verbose = T)

multiomicCCA::compute_enrichment_scores(combined_g, rna_frnn$d_g, atac_frnn$d_g, res$idx)
table(myeloid2@meta.data$celltype[res$idx])

####################

.plot_umap <- function(layout, full_col_vec, idx, main){
  tmp_col_vec <- full_col_vec; tmp_col_vec[-idx] <- "gray"
  size_vec <- rep(0.5, n); size_vec[idx] <- 2
  plot(layout[-idx,1], layout[-idx,2], 
       xlim = range(layout[,1]), ylim = range(layout[,2]),
       asp = T, col = tmp_col_vec[-idx], 
       cex = size_vec[-idx], pch = 16, xlab = "UMAP 1", ylab = "UMAP 2",
       main = main)
  points(layout[idx,1], layout[idx,2],
         col = tmp_col_vec[idx], cex = size_vec[idx])
  invisible()
}

full_col_vec <- scales::hue_pal()(4)[as.numeric(as.factor(myeloid2@meta.data$celltype))]
n <- length(full_col_vec)

png("../../../../out/figures/Writeup14b/Writeup14b_mouseicb_distinct1.png", 
    width = 3000, height = 2000, units = "px", res = 300)
par(mfcol = c(2,3), mar = c(4,4,4,0.1))
layout <- myeloid2[["combined"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, "Common combined,\nInitial")
.plot_umap(layout, full_col_vec, res$idx, "Common combined,\nFinal")

layout <- myeloid2[["distinct"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, "RNA distinct,\nInitial")
.plot_umap(layout, full_col_vec, res$idx, "RNA distinct,\nFinal")

layout <- myeloid2[["distinct2"]]@cell.embeddings
.plot_umap(layout, full_col_vec, start_idx, "ATAC distinct,\nInitial")
.plot_umap(layout, full_col_vec, res$idx, "ATAC distinct,\nFinal")
graphics.off()

