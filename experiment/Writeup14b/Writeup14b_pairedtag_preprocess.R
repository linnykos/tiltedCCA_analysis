rm(list=ls())
library(tidyverse)
library(progressr)
# progressr::handlers(global = TRUE)
# 
# dat <- readr::read_tsv("../../../../data/Pairedtag_mousebrain_RNA-Histone/exprMatrix.tsv")
# 
# # from https://cran.r-project.org/web/packages/progressr/vignettes/progressr-intro.html
# my_fcn <- function(dat2) {
#   p <- progressr::progressor(along = 1:ncol(dat2))
#   y <- purrr::map(dat2, function(x){
#     p()
#     Matrix::Matrix(x, sparse = T)
#   })
# }
# 
# # dat2 is a list of sparse matrix. we will overwrite elements based on one element of said list
# dat2 <- my_fcn(dat[,-1])
# 
# dat3 <- dat2[[1]]
# vec_x <- unlist(pbapply::pblapply(dat2, function(x){x@x}))
# vec_i <- unlist(pbapply::pblapply(dat2, function(x){x@i}))
# vec_p <- c(0, cumsum(pbapply::pbsapply(dat2, function(x){x@p[2]})))
# names(vec_x) <- NULL
# names(vec_i) <- NULL
# names(vec_p) <- NULL
# vec_p <- as.integer(vec_p)
# dat3@x <- vec_x; dat3@i <- vec_i; dat3@p <- vec_p; dat3@Dim[2] <- as.integer(length(dat2))
# dat3@Dimnames[[1]] <- sapply(dat[,1], as.character)
# dat3@Dimnames[[2]] <- colnames(dat)[-1]
# 
# # dat3[100:105,500:505]
# # dat[100:105,501:506]
# 
# dat_tmp <- dat
# dat <- dat3
# 
# save(dat, file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_RNA.RData")

#####################################

dat_1 <- Matrix::readMM("../../../../data/Pairedtag_mousebrain_RNA-Histone/02.Paired-Tag_H3K4me1_DNA_filtered_matrix/matrix.mtx")
tmp <- read.table(file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/02.Paired-Tag_H3K4me1_DNA_filtered_matrix/barcodes.tsv",
                  sep = '\t', header = F)
colname_vec <- tmp[,1]
tmp <- read.table(file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/02.Paired-Tag_H3K4me1_DNA_filtered_matrix/bins.tsv",
                  sep = '\t', header = F)
rowname_vec <- tmp[,1]
dim(dat_1)
class(dat_1)
dat_1 <- as(dat_1, "dgCMatrix")
rownames(dat_1) <- rowname_vec; colnames(dat_1) <- colname_vec
dat_1[1:5,1:5]
all(colnames(dat_1) %in% colnames(dat))

save(dat_1, file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_H3K4me1.RData")


##########

dat_2 <- Matrix::readMM("../../../../data/Pairedtag_mousebrain_RNA-Histone/06.Paired-Tag_H3K9me3_DNA_filtered_matrix/matrix.mtx")
tmp <- read.table(file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/06.Paired-Tag_H3K9me3_DNA_filtered_matrix/barcodes.tsv",
                  sep = '\t', header = F)
colname_vec <- tmp[,1]
tmp <- read.table(file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/06.Paired-Tag_H3K9me3_DNA_filtered_matrix/bins.tsv",
                  sep = '\t', header = F)
rowname_vec <- tmp[,1]
dim(dat_2)
class(dat_2)
dat_2 <- as(dat_2, "dgCMatrix")
rownames(dat_2) <- rowname_vec; colnames(dat_2) <- colname_vec
dat_2[1:5,1:5]
all(colnames(dat_2) %in% colnames(dat))

save(dat_2, file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_H3K9me3.RData")

#############################

dat_rna <- Matrix::readMM("../../../../data/Pairedtag_mousebrain_RNA-Histone/01.Paired-Tag_seq_RNA_filtered_matrix/matrix.mtx")
tmp <- read.table(file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/01.Paired-Tag_seq_RNA_filtered_matrix/barcodes.tsv",
                  sep = '\t', header = F)
colname_vec <- tmp[,1]
tmp <- read.table(file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/01.Paired-Tag_seq_RNA_filtered_matrix/genes.tsv",
                  sep = '\t', header = F)
rowname_vec <- tmp[,1]
dim(dat_rna)
class(dat_rna)
dat_rna <- as(dat_rna, "dgCMatrix")
rownames(dat_rna) <- rowname_vec; colnames(dat_rna) <- colname_vec
dat_rna[1:5,1:5]

save(dat_rna, file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_RNA.RData")

#######################################

# let's make the seurat objects, first for H3K9me3
rm(list=ls())
load("../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_RNA.RData")
load("../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_H3K4me1.RData")
meta_df <- readr::read_tsv("../../../../data/Pairedtag_mousebrain_RNA-Histone/meta.tsv")
meta_df <- as.data.frame(meta_df)
cell_idx <- sapply(colnames(dat_1), function(x){which(colnames(dat_rna) == x)[1]})
cell_idx2 <- sapply(colnames(dat_1), function(x){which(meta_df$Cell_ID == x)[1]})
length(which(is.na(cell_idx)))

dat_rna2 <- dat_rna[,cell_idx]
all(colnames(dat_rna2) == colnames(dat_1))

library(Seurat)
pairedtag <- Seurat::CreateSeuratObject(counts = dat_rna2)
Seurat::DefaultAssay(pairedtag) <- "RNA"
pairedtag[["celltype"]] <- meta_df$Annotation[cell_idx2]
set.seed(10)
pairedtag <- Seurat::SCTransform(pairedtag, verbose = T)
pairedtag <- Seurat::FindVariableFeatures(pairedtag)
pairedtag <- Seurat::RunPCA(pairedtag, verbose = FALSE)
set.seed(10)
pairedtag <- Seurat::RunUMAP(pairedtag, dims = 1:50)

####################

pairedtag[["DNA"]] <- Seurat::CreateAssayObject(counts = dat_1)
Seurat::DefaultAssay(pairedtag) <- "DNA"
pairedtag <- Signac::RunTFIDF(pairedtag)
pairedtag <-  Signac::FindTopFeatures(pairedtag, min.cutoff="q10")
pairedtag <-  Signac::RunSVD(pairedtag)  # question: is this svd run with only the top variable features?
set.seed(10)
pairedtag <- Seurat::RunUMAP(pairedtag, reduction="lsi", dims=1:50,reduction.name="umap.dna", 
                          reduction.key="dnaUMAP_")

save(pairedtag, file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_H3K4me1.RData")

######################

Seurat::DefaultAssay(pairedtag) <- "SCT"
zz1 <- pairedtag[["SCT"]]@data[Seurat::VariableFeatures(pairedtag),]
zz1 <- Matrix::t(zz1)
mean_vec <- sparseMatrixStats::colMeans2(zz1)
sd_vec <- sparseMatrixStats::colSds(zz1)
tmp <- irlba::irlba(zz1, nv = 30, center = mean_vec, scale = sd_vec)
zz1 <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d/max(tmp$d))

Seurat::DefaultAssay(pairedtag) <- "DNA"
zz2 <- pairedtag[["DNA"]]@data[Seurat::VariableFeatures(pairedtag),]
zz2 <- Matrix::t(zz2)
mean_vec <- sparseMatrixStats::colMeans2(zz2)
sd_vec <- sparseMatrixStats::colSds(zz2)
tmp <- irlba::irlba(zz2, nv = 50, center = mean_vec, scale = sd_vec)
zz2 <- multiomicCCA:::.mult_mat_vec(tmp$u[,-1], tmp$d[-1]/max(tmp$d))

rownames(zz1) <- colnames(pairedtag)
rownames(zz2) <- colnames(pairedtag)
zz <- cbind(zz1, zz2)
set.seed(10)
zz <- Seurat::RunUMAP(zz)
pairedtag[["combined"]] <- Seurat::CreateDimReducObject(embedding = zz@cell.embeddings, key = "UMAP", assay = "RNA")

set.seed(10)
zz1_umap <- Seurat::RunUMAP(zz1)
pairedtag[["rnaumap"]] <- Seurat::CreateDimReducObject(embedding = zz1_umap@cell.embeddings, key = "UMAP", assay = "RNA")
set.seed(10)
zz2_umap <- Seurat::RunUMAP(zz2)
pairedtag[["dnaumap"]] <- Seurat::CreateDimReducObject(embedding = zz2_umap@cell.embeddings, key = "UMAP", assay = "RNA")


png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K4me1_base_rna_umap.png", 
    height = 1500, width = 1800, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pairedtag, reduction = "umap", group.by = "celltype", 
                         label = TRUE, repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse Frontal Cortex+Hippo (RNA)\nVanilla UMAP")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1
graphics.off()

png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K4me1_base_histone_umap.png", 
    height = 1500, width = 1800, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pairedtag, reduction = "umap.dna", group.by = "celltype", 
                         label = TRUE, repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse Frontal Cortex+Hippo (Histone)\nVanilla UMAP")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1
graphics.off()


png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K4me1_base_rna_umap2.png", 
    height = 1500, width = 1800, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pairedtag, reduction = "rnaumap", group.by = "celltype", 
                         label = TRUE, repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse Frontal Cortex+Hippo (RNA)\nVanilla UMAP")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1
graphics.off()

png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K4me1_base_histone_umap2.png", 
    height = 1500, width = 1800, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pairedtag, reduction = "dnaumap", group.by = "celltype", 
                         label = TRUE, repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse Frontal Cortex+Hippo (Histone)\nVanilla UMAP")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1
graphics.off()

png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K4me1_base_both_umap.png", 
    height = 1500, width = 1800, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pairedtag, reduction = "combined", group.by = "celltype", 
                         label = TRUE, repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse Frontal Cortex+Hippo (RNA+Histone)\nVanilla UMAP")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1
graphics.off()

#####################################################3

# let's make the seurat objects, now for H3K9me3
rm(list=ls())
load("../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_RNA.RData")
load("../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_H3K9me3.RData")
meta_df <- readr::read_tsv("../../../../data/Pairedtag_mousebrain_RNA-Histone/meta.tsv")
meta_df <- as.data.frame(meta_df)
cell_idx <- sapply(colnames(dat_2), function(x){which(colnames(dat_rna) == x)[1]})
cell_idx2 <- sapply(colnames(dat_2), function(x){which(meta_df$Cell_ID == x)[1]})
length(which(is.na(cell_idx)))

dat_rna2 <- dat_rna[,cell_idx]
all(colnames(dat_rna2) == colnames(dat_2))

library(Seurat)
pairedtag <- Seurat::CreateSeuratObject(counts = dat_rna2)
Seurat::DefaultAssay(pairedtag) <- "RNA"
pairedtag[["celltype"]] <- meta_df$Annotation[cell_idx2]
set.seed(10)
pairedtag <- Seurat::SCTransform(pairedtag, verbose = T)
pairedtag <- Seurat::FindVariableFeatures(pairedtag)
pairedtag <- Seurat::RunPCA(pairedtag, verbose = FALSE)
set.seed(10)
pairedtag <- Seurat::RunUMAP(pairedtag, dims = 1:50)

####################

pairedtag[["DNA"]] <- Seurat::CreateAssayObject(counts = dat_2)
Seurat::DefaultAssay(pairedtag) <- "DNA"
pairedtag <- Signac::RunTFIDF(pairedtag)
pairedtag <-  Signac::FindTopFeatures(pairedtag, min.cutoff="q10")
pairedtag <-  Signac::RunSVD(pairedtag)  # question: is this svd run with only the top variable features?
set.seed(10)
pairedtag <- Seurat::RunUMAP(pairedtag, reduction="lsi", dims=1:50,reduction.name="umap.dna", 
                             reduction.key="dnaUMAP_")

#############

Seurat::DefaultAssay(pairedtag) <- "SCT"
zz1 <- pairedtag[["SCT"]]@data[Seurat::VariableFeatures(pairedtag),]
zz1 <- Matrix::t(zz1)
mean_vec <- sparseMatrixStats::colMeans2(zz1)
sd_vec <- sparseMatrixStats::colSds(zz1)
tmp <- irlba::irlba(zz1, nv = 50)#, center = mean_vec) #, scale = sd_vec)
zz1 <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d/max(tmp$d))

Seurat::DefaultAssay(pairedtag) <- "DNA"
zz2 <- pairedtag[["DNA"]]@data[Seurat::VariableFeatures(pairedtag),]
zz2 <- Matrix::t(zz2)
mean_vec <- sparseMatrixStats::colMeans2(zz2)
sd_vec <- sparseMatrixStats::colSds(zz2)
tmp <- irlba::irlba(zz2, nv = 50)# , center = mean_vec) #, scale = sd_vec)
zz2 <- multiomicCCA:::.mult_mat_vec(tmp$u, tmp$d/max(tmp$d))

rownames(zz1) <- colnames(pairedtag)
rownames(zz2) <- colnames(pairedtag)
zz <- cbind(zz1, zz2)
set.seed(10)
zz <- Seurat::RunUMAP(zz)
pairedtag[["combined"]] <- Seurat::CreateDimReducObject(embedding = zz@cell.embeddings, key = "UMAP", assay = "RNA")

set.seed(10)
zz1_umap <- Seurat::RunUMAP(zz1)
pairedtag[["rnaumap"]] <- Seurat::CreateDimReducObject(embedding = zz1_umap@cell.embeddings, key = "UMAP", assay = "RNA")
set.seed(10)
zz2_umap <- Seurat::RunUMAP(zz2)
pairedtag[["dnaumap"]] <- Seurat::CreateDimReducObject(embedding = zz2_umap@cell.embeddings, key = "UMAP", assay = "RNA")


png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K9me3_base_rna_umap.png", 
    height = 1500, width = 1800, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pairedtag, reduction = "umap", group.by = "celltype", 
                         label = TRUE, repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse Frontal Cortex+Hippo (RNA)\nH3K9me3: Vanilla UMAP")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1
graphics.off()

png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K9me3_base_histone_umap.png", 
    height = 1500, width = 1800, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pairedtag, reduction = "umap.dna", group.by = "celltype", 
                         label = TRUE, repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse Frontal Cortex+Hippo (Histone)\nH3K9me3: Vanilla UMAP")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1
graphics.off()


png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K9me3_base_rna_umap2.png", 
    height = 1500, width = 1800, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pairedtag, reduction = "rnaumap", group.by = "celltype", 
                         label = TRUE, repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse Frontal Cortex+Hippo (RNA)\nH3K9me3: Vanilla UMAP2")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1
graphics.off()

png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K9me3_base_histone_umap2.png", 
    height = 1500, width = 1800, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pairedtag, reduction = "dnaumap", group.by = "celltype", 
                         label = TRUE, repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse Frontal Cortex+Hippo (Histone)\nH3K9me3: Vanilla UMAP2")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1
graphics.off()

png("../../../../out/figures/Writeup14b/Writeup14b_pairedtag_mousecortex_H3K9me3_base_both_umap.png", 
    height = 1500, width = 1800, units = "px", res = 300)
plot1 <- Seurat::DimPlot(pairedtag, reduction = "combined", group.by = "celltype", 
                         label = TRUE, repel = TRUE, label.size = 2.5) 
plot1 <- plot1 + ggplot2::ggtitle("Mouse Frontal Cortex+Hippo (RNA+Histone)\nH3K9me3: Vanilla UMAP")
plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
plot1
graphics.off()

