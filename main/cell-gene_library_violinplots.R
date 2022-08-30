rm(list=ls())

library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/citeseq_bm25_preprocessed.RData")
bm25_celldepth <- bm$nCount_RNA
bm25_celldepth_ab <- bm$nCount_ADT
bm25_genedepth <- Matrix::rowSums(bm[["RNA"]]@counts[bm[["RNA"]]@var.features,])
bm25_genedepth <- bm25_genedepth/ncol(bm)
bm25_abdepth <- Matrix::rowSums(bm[["ADT"]]@counts[bm[["ADT"]]@var.features,])
bm25_abdepth <- bm25_abdepth/ncol(bm)
round(quantile(bm25_celldepth),2)
round(quantile(bm25_celldepth_ab),2)
round(quantile(bm25_genedepth),2)
round(quantile(bm25_abdepth),2)

load("../../../out/main/abseq_bm97_preprocessed.RData")
bm97_celldepth <- bm$nCount_RNA
bm97_celldepth_ab <- bm$nCount_AB
bm97_genedepth <- Matrix::rowSums(bm[["RNA"]]@counts[bm[["RNA"]]@var.features,])
bm97_genedepth <- bm97_genedepth/ncol(bm)
bm97_abdepth <- Matrix::rowSums(bm[["AB"]]@counts[bm[["AB"]]@var.features,])
bm97_abdepth <- bm97_abdepth/ncol(bm)
round(quantile(bm97_celldepth),2)
round(quantile(bm97_celldepth_ab),2)
round(quantile(bm97_genedepth),2)
round(quantile(bm97_abdepth),2)

load("../../../out/main/abseq_bm97Ref_preprocessed.RData")
bm97ref_celldepth <- bm$nCount_RNA
bm97ref_celldepth_ab <- bm$nCount_AB
bm97ref_genedepth <- Matrix::rowSums(bm[["RNA"]]@counts[bm[["RNA"]]@var.features,])
bm97ref_genedepth <- bm97ref_genedepth/ncol(bm)
bm97ref_abdepth <- Matrix::rowSums(bm[["AB"]]@counts[bm[["AB"]]@var.features,])
bm97ref_abdepth <- bm97ref_abdepth/ncol(bm)
round(quantile(bm97ref_celldepth),2)
round(quantile(bm97ref_celldepth_ab),2)
round(quantile(bm97ref_genedepth),2)
round(quantile(bm97ref_abdepth),2)

load("../../../out/main/citeseq_pbmc224_preprocessed.RData")
pbmc224_celldepth <- pbmc$nCount_RNA
pbmc224_celldepth_ab <- pbmc$nCount_ADT
pbmc224_genedepth <- Matrix::rowSums(pbmc[["SCT"]]@counts[pbmc[["SCT"]]@var.features,])
pbmc224_genedepth <- pbmc224_genedepth/ncol(pbmc)
pbmc224_abdepth <- Matrix::rowSums(pbmc[["ADT"]]@counts)
pbmc224_abdepth <- pbmc224_abdepth/ncol(pbmc)
round(quantile(pbmc224_celldepth),2)
round(quantile(pbmc224_celldepth_ab),2)
round(quantile(pbmc224_genedepth),2)
round(quantile(pbmc224_abdepth),2)

load("../../../out/main/10x_pbmc_preprocessed.RData")
pbmc10x_celldepth <- pbmc$nCount_RNA
pbmc10x_celldepth_atac <- pbmc$nCount_ATAC
pbmc10x_genedepth <- Matrix::rowSums(pbmc[["RNA"]]@counts[pbmc[["SCT"]]@var.features,])
pbmc10x_genedepth <- pbmc10x_genedepth/ncol(pbmc)
round(quantile(pbmc10x_celldepth),2)
round(quantile(pbmc10x_celldepth_atac),2)
round(quantile(pbmc10x_genedepth),2)

load("../../../out/main/10x_greenleaf_preprocessed.RData")
greenleaf_celldepth <- greenleaf$nCount_RNA
greenleaf_celldepth_atac <- greenleaf$nCount_ATAC
greenleaf_genedepth <- Matrix::rowSums(greenleaf[["RNA"]]@counts[greenleaf[["SCT"]]@var.features,])
greenleaf_genedepth <- greenleaf_genedepth/ncol(greenleaf)
round(quantile(greenleaf_celldepth),2)
round(quantile(greenleaf_celldepth_atac),2)
round(quantile(greenleaf_genedepth),2)

load("../../../out/main/10x_mouseembryo_preprocessed.RData")
mbrain_celldepth <- mbrain$nCount_RNA
mbrain_celldepth_atac <- mbrain$nCount_ATAC
mbrain_genedepth <- Matrix::rowSums(mbrain[["RNA"]]@counts[mbrain[["SCT"]]@var.features,])
mbrain_genedepth <- mbrain_genedepth/ncol(mbrain)
round(quantile(mbrain_celldepth),2)
round(quantile(mbrain_celldepth_atac),2)
round(quantile(mbrain_genedepth),2)

load("../../../out/main/10x_reik_preprocessed.RData")
reik_celldepth <- reik$nCount_RNA
reik_celldepth_atac <- reik$nCount_ATAC
reik_genedepth <- Matrix::rowSums(reik[["RNA"]]@counts[reik[["RNA"]]@var.features,])
reik_genedepth <- reik_genedepth/ncol(reik)
round(quantile(reik_celldepth),2)
round(quantile(reik_celldepth_atac),2)
round(quantile(reik_genedepth),2)

################

celltype_depth <- data.frame(
  len = c(bm25_celldepth, bm97_celldepth, bm97ref_celldepth, pbmc224_celldepth,
          pbmc10x_celldepth, greenleaf_celldepth, mbrain_celldepth, reik_celldepth),
  dataset = c(rep("BM25", length(bm25_celldepth)), rep("BM97", length(bm97_celldepth)),
              rep("BM97Ref", length(bm97ref_celldepth)), rep("PBMC224", length(pbmc224_celldepth)),
              rep("PBMC10x", length(pbmc10x_celldepth)), rep("Greenleaf", length(greenleaf_celldepth)),
              rep("Mbrain", length(mbrain_celldepth)), rep("Reik", length(reik_celldepth)))
)
celltype_depth$len <- log10(celltype_depth$len+1)

# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
plot1 <- ggplot2::ggplot(celltype_depth, ggplot2::aes(x=dataset, y=len)) + 
  ggplot2::geom_violin(trim=FALSE, fill="gray")+
  ggplot2::geom_boxplot(width=0.1, outlier.shape = NA)+
  ggplot2::theme_classic()
ggplot2::ggsave(filename = paste0("../../../out/figures/main/celldepth_violin.png"),
                plot1, device = "png", width = 10, height = 3, units = "in",
                dpi = 500)

################

celltype_depth_ab <- data.frame(
  len = c(bm25_celldepth_ab, bm97_celldepth_ab, bm97ref_celldepth_ab, pbmc224_celldepth_ab),
  dataset = c(rep("BM25", length(bm25_celldepth_ab)), rep("BM97", length(bm97_celldepth_ab)),
              rep("BM97Ref", length(bm97ref_celldepth_ab)), rep("PBMC224", length(pbmc224_celldepth_ab)))
)
celltype_depth_ab$len <- log10(celltype_depth_ab$len+1)

plot1 <- ggplot2::ggplot(celltype_depth_ab, ggplot2::aes(x=dataset, y=len)) + 
  ggplot2::geom_violin(trim=FALSE, fill="gray")+
  ggplot2::geom_boxplot(width=0.1, outlier.shape = NA)+
  ggplot2::theme_classic()
ggplot2::ggsave(filename = paste0("../../../out/figures/main/celldepth_ab_violin.png"),
                plot1, device = "png", width = 5, height = 3, units = "in",
                dpi = 500)

################

celltype_depth_atac <- data.frame(
  len = c(pbmc10x_celldepth_atac, greenleaf_celldepth_atac, mbrain_celldepth_atac, reik_celldepth_atac),
  dataset = c(rep("PBMC10x", length(pbmc10x_celldepth_atac)), rep("Greenleaf", length(greenleaf_celldepth_atac)),
              rep("Mbrain", length(mbrain_celldepth_atac)), rep("Reik", length(reik_celldepth_atac)))
)
celltype_depth_atac$len <- log10(celltype_depth_atac$len+1)

plot1 <- ggplot2::ggplot(celltype_depth_atac, ggplot2::aes(x=dataset, y=len)) + 
  ggplot2::geom_violin(trim=FALSE, fill="gray")+
  ggplot2::geom_boxplot(width=0.1, outlier.shape = NA)+
  ggplot2::theme_classic()
ggplot2::ggsave(filename = paste0("../../../out/figures/main/celldepth_atac_violin.png"),
                plot1, device = "png", width = 5, height = 3, units = "in",
                dpi = 500)

################

gene_depth <- data.frame(
  len = c(bm25_genedepth, bm97_genedepth, bm97ref_genedepth, pbmc224_genedepth,
          pbmc10x_genedepth, greenleaf_genedepth, mbrain_genedepth, reik_genedepth),
  dataset = c(rep("BM25", length(bm25_genedepth)), rep("BM97", length(bm97_genedepth)),
              rep("BM97Ref", length(bm97ref_genedepth)), rep("PBMC224", length(pbmc224_genedepth)),
              rep("PBMC10x", length(pbmc10x_genedepth)), rep("Greenleaf", length(greenleaf_genedepth)),
              rep("Mbrain", length(mbrain_genedepth)), rep("Reik", length(reik_genedepth)))
)
gene_depth$len <- log10(gene_depth$len+1)

# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
plot1 <- ggplot2::ggplot(gene_depth, ggplot2::aes(x=dataset, y=len)) + 
  ggplot2::geom_violin(trim=FALSE, fill="gray")+
  ggplot2::geom_boxplot(width=0.1, outlier.shape = NA)+
  ggplot2::theme_classic()
ggplot2::ggsave(filename = paste0("../../../out/figures/main/genedepth_violin.png"),
                plot1, device = "png", width = 10, height = 3, units = "in",
                dpi = 500)


################

ab_depth <- data.frame(
  len = c(bm25_abdepth, bm97_abdepth, bm97ref_abdepth, pbmc224_abdepth),
  dataset = c(rep("BM25", length(bm25_abdepth)), rep("BM97", length(bm97_abdepth)),
              rep("BM97Ref", length(bm97ref_abdepth)), rep("PBMC224", length(pbmc224_abdepth)))
)
ab_depth$len <- log10(ab_depth$len+1)

# http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
plot1 <- ggplot2::ggplot(ab_depth, ggplot2::aes(x=dataset, y=len)) + 
  ggplot2::geom_violin(trim=FALSE, fill="gray")+
  ggplot2::geom_boxplot(width=0.1, outlier.shape = NA)+
  ggplot2::theme_classic()
ggplot2::ggsave(filename = paste0("../../../out/figures/main/abdepth_violin.png"),
                plot1, device = "png", width = 5, height = 3, units = "in",
                dpi = 500)


################

date_of_run <- Sys.time()
session_info <- devtools::session_info()
set.seed(10)

save(gene_depth, ab_depth,
     celltype_depth, celltype_depth_ab, celltype_depth_atac, 
     date_of_run, session_info,
     file = "../../../out/main/cell-gene_library.RData")

######################

sapply(sort(unique(celltype_depth[,2])), function(x){
  10^quantile(celltype_depth[which(celltype_depth[,2] == x),1])-1
})
