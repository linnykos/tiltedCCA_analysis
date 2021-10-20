rm(list=ls())
genotype_est <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_est.rds")
genotype_values <- readRDS("../../../../data/Chiyun_SNU601/SNU601/genotype_values_h1_h2_c20.rds")
seg_table_filtered <- readRDS("../../../../data/Chiyun_SNU601/SNU601/seg_table_filtered.rds")
mat <- Matrix::readMM("../../../../data/Chiyun_SNU601/SNU601/matrix.mtx")
features <- read.csv("../../../../data/Chiyun_SNU601/SNU601/features.txt", header = F)
barcodes <- read.csv("../../../../data/Chiyun_SNU601/SNU601/barcodes.txt", header = F)
clone_assign <- readRDS("../../../../data/Chiyun_SNU601/SNU601/cloneAssign.rds")
qc_per_barcode <- read.csv("../../../../data/Chiyun_SNU601/SNU601/SNU601.MACS2.qc_per_barcode.txt", 
                           header = T, sep = "\t")

rownames(mat) <- features[,1]
colnames(mat) <- barcodes[,1]
mat <- mat[,names(clone_assign)]

##############

all(colnames(mat) %in% qc_per_barcode$bc)
idx_vec <- sapply(colnames(mat), function(barcode){
  which(qc_per_barcode$bc == barcode)[1]
})
qc_per_barcode <- qc_per_barcode[idx_vec,]
rownames(qc_per_barcode) <- qc_per_barcode[,1]
qc_per_barcode <- qc_per_barcode[,-1]

library(Seurat); library(Signac); library("EnsDb.Hsapiens.v75")

chrom_assay <- Signac::CreateChromatinAssay(counts = mat,
                                            sep = c("-", "-"),
                                            genome = "hg19",
                                            fragments = '../../../../data/Chiyun_SNU601/SNU601/fragments.tsv.gz',
                                            min.cells = 1)
SNU <- Seurat::CreateSeuratObject(
  counts = chrom_assay,
  assay = "atac",
  meta.data = qc_per_barcode
)
SNU[["clone"]] <- clone_assign

###################
# from https://satijalab.org/signac/articles/pbmc_vignette.html
# extract gene annotations from EnsDb
annotations <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
# change to UCSC style since the data was mapped to hg19
GenomeInfoDb::seqlevelsStyle(annotations) <- 'UCSC'
GenomeInfoDb::genome(annotations) <- "hg19"
# add the gene information to the object
Signac::Annotation(SNU) <- annotations

#############

SNU[["atac"]]
GenomicRanges::granges(SNU)
SNU <- Signac::NucleosomeSignal(object = SNU)
SNU <- Signac::TSSEnrichment(object = SNU, fast = FALSE)

############

SNU <- Signac::FindTopFeatures(SNU, min.cutoff = 'q10')
mat <- mat[SNU[["atac"]]@var.features,]
set.seed(10)
svd_res <- irlba::irlba(mat, nv = 200, verbose = T)
v_mat <- multiomicCCA:::.mult_mat_vec(svd_res$v, svd_res$d)
for(j in 1:ncol(v_mat)){
  df <- data.frame(y = v_mat[,j], library_size = SNU@meta.data[,"nCount_atac"])
  lm_res <- stats::lm(y ~ ., data = df)
  v_mat[,j] <- stats::residuals(lm_res)
}
mat2 <- tcrossprod(svd_res$u, v_mat)
colnames(mat2) <- colnames(mat)
rownames(mat2) <- rownames(mat)
SNU[["atac"]]@data <- mat2

SNU[["lsi"]] <- Signac::RunSVD(SNU[["atac"]]@data,
                               assay = "atac")
set.seed(10)
SNU <- Seurat::RunUMAP(SNU, reduction = 'lsi', dims = 1:50, 
                       reduction.name = "umap.atac", 
                       reduction.key = "atacUMAP_")
rownames(SNU[["umap.atac"]]@cell.embeddings) <- rownames(SNU@meta.data)
save(SNU, file = "../../../../out/Writeup14e/Writeup14e_SNU_atac_exploratory.RData")

#####################


group_vec <- c("clone")
for(group in group_vec){
  plot1 <- Seurat::DimPlot(SNU, reduction = "umap.atac",
                           group.by = group, label = TRUE,
                           repel = TRUE, label.size = 2.5)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (ATAC)\n",group))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_atac_exploration_", group, ".png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
}

group_vec <- c("nCount_atac", "nFeature_atac", "total_frags",
               "frac_peak", "frac_promoter", "frac_tss",
               "frac_enhancer", "nucleosome_signal",
               "nucleosome_percentile", "TSS.enrichment", "TSS.percentile")
for(group in group_vec){
  plot1 <- Seurat::FeaturePlot(SNU, reduction = "umap.atac",
                               features = group)
  plot1 <- plot1 + ggplot2::ggtitle(paste0("SNU (ATAC)\n",group))
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../../out/figures/Writeup14e/Writeup14e_SNU_atac_exploration_", group, ".png"),
                  plot1, device = "png", width = 5, height = 5, units = "in")
  
}


