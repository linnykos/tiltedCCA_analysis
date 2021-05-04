rm(list=ls())
load("../../../../data/Pairedtag_mousebrain_RNA-Histone/dat.RData")
library(Matrix)

dat <- t(dat)
sd_vec <- sparseMatrixStats::colSds(dat)
mean_vec <- sparseMatrixStats::colMeans2(dat)

set.seed(10)
svd_res <- irlba::irlba(dat, nv = 50, center = mean_vec, scale = sd_vec)
svd_res$d
diff(svd_res$d)/svd_res$d[-1]
k <- 19
dim_red <- .mult_mat_vec(svd_res$u[,1:k], svd_res$d[1:k])
set.seed(10)
umap_res <- Seurat::RunUMAP(dim_red, verbose = T)

#############

meta_df <- readr::read_tsv("../../../../data/Pairedtag_mousebrain_RNA-Histone/meta.tsv")
meta_df <- as.data.frame(meta_df)
obj <- Seurat::CreateSeuratObject(t(dat), meta.data = meta_df)

obj[["asdf"]] <- umap_res
png("../../../../out/figures/Writeup14/Writeup14_pairedtag_eda.png", height = 1500, width = 1500, units = "px", res = 300)
plot1 <- Seurat::DimPlot(obj, reduction = 'asdf', group.by = 'Annotation', label = TRUE, 
                         repel = TRUE, label.size = 2.5) + Seurat::NoLegend()
plot1 + ggplot2::ggtitle("EDA (Paired-Tag)")
graphics.off()
