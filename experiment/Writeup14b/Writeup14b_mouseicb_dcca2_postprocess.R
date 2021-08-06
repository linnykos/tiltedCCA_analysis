rm(list=ls())
load("../../../../out/Writeup14b/Writeup14b_mouseicb_dcca.RData")
load("../../../../out/Writeup14b/Writeup14b_mouseicb_dcca_tmp.RData")

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res, verbose = T)
mat_1_denoised <- dcca_decomp$common_mat_1 + dcca_decomp$distinct_mat_1
mat_2_denoised <- dcca_decomp$common_mat_2 + dcca_decomp$distinct_mat_2

##########

c_eig <- myeloid2[["clap"]]@cell.embeddings
d_eig <- myeloid2[["dlap"]]@cell.embeddings
e_eig <- myeloid2[["elap"]]@cell.embeddings

idx <- which(rownames(dcca_res$svd_1$v) == "Stat1")
c_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,idx], c_eig)
d_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,idx], d_eig)
e_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,idx], e_eig)

tmp <- multiomicCCA::plot_laplacian(myeloid2, var_name = colnames(mat_1_denoised)[idx],
                             prefix = "RNA", e_vec = mat_1_denoised[,idx],
                             c_vec = dcca_decomp$common_mat_1[,idx],
                             d_vec = dcca_decomp$distinct_mat_1[,idx],
                             e_res = e_res, c_res = c_res, d_res = d_res,
                             reduction_1 = "everything", reduction_2 = "combined", reduction_3 = "distinct")
tmp2 <- cowplot::plot_grid(tmp[[1]], tmp[[2]], tmp[[3]], tmp[[4]], tmp[[5]], tmp[[6]])
cowplot::save_plot(filename = "../../../../out/figures/Writeup14b/Writeup14b_mouseicb_dcca_rna_stat1_umap.png",
                   tmp2, ncol = 3, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")


val_vec <- sapply(gene_smoothed, function(x){
  (x$d_variance - x$c_variance)/x$e_variance
})
name_vec <- colnames(mat_1_denoised)
p1 <- length(name_vec)
factor_vec <- rep(0, p1); 
factor_vec[which(sapply(gene_smoothed, function(x){min(x$c_r2, x$d_r2) > 0.6}))] <- 2
factor_vec[idx] <- 1; factor_vec <- as.factor(factor_vec)
col_vec <- c("black", "red", "green"); names(col_vec) <- c("0", "1", "2")
tmp <- plot_laplacian_variables(val_vec, name_vec, factor_vec, col_vec,
                                              ylab = "Distinct-common (normalized)", main = "RNA enrichment",
                                text_cex = 4)
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14b/Writeup14b_mouseicb_dcca_rna_var_enrichment.png",
                tmp, device = "png", width = 3, height = 3, units = "in")

idx1 <- which(factor_vec == 2)
idx2 <- which.max(val_vec[idx1])
idx <- idx1[idx2]
c_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,idx], c_eig)
d_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,idx], d_eig)
e_res <- multiomicCCA::compute_smooth_signal(mat_1_denoised[,idx], e_eig)

tmp <- multiomicCCA::plot_laplacian(myeloid2, var_name = colnames(mat_1_denoised)[idx],
                                    prefix = "RNA", e_vec = mat_1_denoised[,idx],
                                    c_vec = dcca_decomp$common_mat_1[,idx],
                                    d_vec = dcca_decomp$distinct_mat_1[,idx],
                                    e_res = e_res, c_res = c_res, d_res = d_res,
                                    reduction_1 = "everything", reduction_2 = "combined", reduction_3 = "distinct")
tmp2 <- cowplot::plot_grid(tmp[[1]], tmp[[2]], tmp[[3]], tmp[[4]], tmp[[5]], tmp[[6]])
cowplot::save_plot(filename = "../../../../out/figures/Writeup14b/Writeup14b_mouseicb_dcca_rna_highest_umap.png",
                   tmp2, ncol = 3, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")

###########################
###########################

myeloid <- readRDS("../../../../data/ICB_mouse/Object_seurat_2000_each.rds")

getRegion <- function(vicinity.size, TSS){
  upstream.size <- vicinity.size[1]
  downstream.size <- vicinity.size[2]
  if (upstream.size == 0) {
    upstream <- TSS
  } else if (upstream.size > 0) {
    upstream <- GenomicRanges::flank(TSS, width=upstream.size, start=T) # upstream excluding the TSS
  }
  if (downstream.size == 0) {
    downstream <- TSS
  } else if (downstream.size > 0) {
    downstream <- GenomicRanges::flank(TSS, width=downstream.size, start=F) # downstream excluding the TSS
  }
  vicinity <- GenomicRanges::punion(upstream, downstream, fill.gap=T)
  BiocGenerics::start(vicinity) <- pmax(0, BiocGenerics::start(vicinity))
  vicinity$gene_name <- TSS$gene_name
  return(vicinity)
}

peaks_gr <- GenomicRanges::granges(myeloid@assays$ATAC)
annotations <- myeloid@assays$ATAC@annotation

temp <- Signac::LookupGeneCoords(myeloid, "Stat1")
if(is.null(temp)) return(numeric(0))
GenomeInfoDb::seqlevelsStyle(temp)="UCSC"  # just so that the seqlevelsStyle matches up with peaks_gr
temp2 <- annotations[which(annotations$gene_name=="Stat1" & annotations$type=="cds")] #cds stands for "coding sequence"
if(length(temp2)==0) return(numeric(0))
gene_strand <- BiocGenerics::strand(temp2)[1]
TSS_position <- ifelse(gene_strand == "+", BiocGenerics::start(temp), BiocGenerics::end(temp))
TSS <- GenomicRanges::GRanges(seqnames = GenomeInfoDb::seqnames(temp),
                              ranges = IRanges::IRanges(start = TSS_position, width = 1),
                              strand = BiocGenerics::strand(temp),
                              gene_name = "Stat1")
upstream_gr <- getRegion(c(50000, 0), TSS)
downstream_gr <- getRegion(c(0, 50000), TSS)

# populate the peaks counts array.
peaks_upstream <- GenomicRanges::countOverlaps(peaks_gr, upstream_gr)>0
peaks_downstream <- GenomicRanges::countOverlaps(peaks_gr, downstream_gr)>0

target_idx <- c(sort(unique(c(which(peaks_upstream), which(peaks_downstream)))))
target_names <- rownames(myeloid[["ATAC"]])[target_idx]
rm(list = c("myeloid", "peaks_gr"))

###########################
###########################

zz <- which(sapply(atac_smoothed, length) > 0)
atac_smoothed2 <- atac_smoothed[zz]

val_vec <- sapply(atac_smoothed2, function(x){
  (x$d_variance - x$c_variance)/x$e_variance
})
name_vec <- colnames(mat_2_denoised)[zz]; name_vec2 <- name_vec
p2 <- length(name_vec)
factor_vec <- rep(0, p2)
factor_vec[which(sapply(atac_smoothed2, function(x){min(x$c_r2, x$d_r2) > 0.5}))] <- 2
factor_vec[which(name_vec %in% target_names)] <- 1; factor_vec <- as.factor(factor_vec)
col_vec <- c("black", "red", "green"); names(col_vec) <- c("0", "1", "2")
name_vec2[which(name_vec %in% target_names)[1]] <- "Peaks for Stat1"
name_vec2[which(name_vec %in% target_names)[-1]] <- ""
tmp <- plot_laplacian_variables(val_vec, name_vec2, factor_vec, col_vec,
                                ylab = "Distinct-common (normalized)", main = "ATAC enrichment",
                                text_cex = 3)
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14b/Writeup14b_mouseicb_dcca_atac_var_enrichment.png",
                tmp, device = "png", width = 3, height = 3, units = "in")

idx1 <- which(factor_vec == 2)
idx2 <- which.max(val_vec[idx1])
idx <- idx1[idx2]

c_eig2 <- myeloid2[["clap2"]]@cell.embeddings
d_eig2 <- myeloid2[["dlap2"]]@cell.embeddings
e_eig2 <- myeloid2[["elap2"]]@cell.embeddings

c_res <- multiomicCCA::compute_smooth_signal(mat_2_denoised[,zz[idx]], c_eig2)
d_res <- multiomicCCA::compute_smooth_signal(mat_2_denoised[,zz[idx]], d_eig2)
e_res <- multiomicCCA::compute_smooth_signal(mat_2_denoised[,zz[idx]], e_eig2)

tmp <- multiomicCCA::plot_laplacian(myeloid2, var_name = colnames(mat_2_denoised)[zz[idx]],
                                    prefix = "ATAC", e_vec = mat_2_denoised[,zz[idx]],
                                    c_vec = dcca_decomp$common_mat_2[,zz[idx]],
                                    d_vec = dcca_decomp$distinct_mat_2[,zz[idx]],
                                    e_res = e_res, c_res = c_res, d_res = d_res,
                                    reduction_1 = "everything2", reduction_2 = "combined", reduction_3 = "distinct2")
tmp2 <- cowplot::plot_grid(tmp[[1]], tmp[[2]], tmp[[3]], tmp[[4]], tmp[[5]], tmp[[6]])
cowplot::save_plot(filename = "../../../../out/figures/Writeup14b/Writeup14b_mouseicb_dcca_atac_highest_umap.png",
                   tmp2, ncol = 3, nrow = 2, base_height = 3.5, base_asp = 4/3, device = "png")

