rm(list=ls())
library(Seurat)
library(Matrix)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

load("../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_RNA.RData")
dna_dat <- Matrix::readMM("../../../../data/Pairedtag_mousebrain_RNA-Histone/07.Paired-seq_DNA_filtered_matrix/matrix.mtx")
dim(dna_dat)
dna_dat <- as(dna_dat, "dgCMatrix")
cell_id <- readr::read_tsv("../../../../data/Pairedtag_mousebrain_RNA-Histone/07.Paired-seq_DNA_filtered_matrix/barcodes.tsv",
                           col_names = F)
cell_id <- as.data.frame(cell_id)
length(which(cell_id[,1] %in% colnames(dat)))/nrow(cell_id)

cell_id2 <- readr::read_tsv("../../../../data/Pairedtag_mousebrain_RNA-Histone/02.Paired-Tag_H3K4me1_DNA_filtered_matrix/barcodes.tsv",
                            col_names = F)
cell_id2 <- as.data.frame(cell_id2)
length(which(cell_id2[,1] %in% colnames(dat)))/nrow(cell_id2)
length(which(cell_id2[,1] %in% cell_id[,1]))
