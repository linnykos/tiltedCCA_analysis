rm(list=ls())
library(tidyverse)
library(progressr)
library(Seurat)
library(Signac)

folder_names <- c("02.Paired-Tag_H3K4me1_DNA_filtered_matrix",
                  "03.Paired-Tag_H3K4me3_DNA_filtered_matrix",
                  "04.Paired-Tag_H3K27ac_DNA_filtered_matrix",
                  "05.Paired-Tag_H3K27me3_DNA_filtered_matrix",
                  "06.Paired-Tag_H3K9me3_DNA_filtered_matrix")
histone_names <- c("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K9me3")

for(i in 1:length(folder_names)){
  print(paste0("Formatting histone ", histone_names[i]))
  
  dat_1 <- Matrix::readMM(paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/",
                                 folder_names[i],
                                 "/matrix.mtx"))
  tmp <- read.table(file = paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/",
                                  folder_names[i],
                                  "/barcodes.tsv"),
                    sep = '\t', header = F)
  colname_vec <- tmp[,1]
  tmp <- read.table(file = paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/",
                                  folder_names[i],
                                  "/bins.tsv"),
                    sep = '\t', header = F)
  rowname_vec <- tmp[,1]
  
  dat_1 <- as(dat_1, "dgCMatrix")
  rownames(dat_1) <- rowname_vec; colnames(dat_1) <- colname_vec
  
  save(dat_1, 
       file = paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_",
                     histone_names[i],
                     ".RData"))
}

#############################
print("Formatting RNA")
dat_rna <- Matrix::readMM("../../../../data/Pairedtag_mousebrain_RNA-Histone/01.Paired-Tag_seq_RNA_filtered_matrix/matrix.mtx")
tmp <- read.table(file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/01.Paired-Tag_seq_RNA_filtered_matrix/barcodes.tsv",
                  sep = '\t', header = F)
colname_vec <- tmp[,1]
tmp <- read.table(file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/01.Paired-Tag_seq_RNA_filtered_matrix/genes.tsv",
                  sep = '\t', header = F)
rowname_vec <- tmp[,1]
dat_rna <- as(dat_rna, "dgCMatrix")
rownames(dat_rna) <- rowname_vec; colnames(dat_rna) <- colname_vec

save(dat_rna, file = "../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_RNA.RData")

#######################################

# let's make the seurat objects
ls_vec <- ls()
ls_vec <- ls_vec[!ls_vec %in% c("folder_names", "histone_names")]
rm(list = ls_vec)

for(i in 1:length(histone_names)){
  print(paste0("Working on Seurat-processing for histone ", histone_names[i]))
  
  load("../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_RNA.RData")
  load(paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/preprocessed_",
              histone_names[i], ".RData"))
  
  meta_df <- readr::read_tsv("../../../../data/Pairedtag_mousebrain_RNA-Histone/meta.tsv")
  meta_df <- as.data.frame(meta_df)
  cell_idx <- sapply(colnames(dat_1), function(x){which(colnames(dat_rna) == x)[1]})
  cell_idx2 <- sapply(colnames(dat_1), function(x){which(meta_df$Cell_ID == x)[1]})
  
  dat_rna2 <- dat_rna[,cell_idx]
  
  pairedtag <- Seurat::CreateSeuratObject(counts = dat_rna2)
  Seurat::DefaultAssay(pairedtag) <- "RNA"
  pairedtag[["celltype"]] <- meta_df$Annotation[cell_idx2]
  set.seed(10)
  pairedtag <- Seurat::SCTransform(pairedtag, verbose = T)
  pairedtag <- Seurat::FindVariableFeatures(pairedtag)
  pairedtag <- Seurat::RunPCA(pairedtag, verbose = FALSE)
  set.seed(10)
  pairedtag <- Seurat::RunUMAP(pairedtag, dims = 1:50)
  
  
  pairedtag[["DNA"]] <- Seurat::CreateAssayObject(counts = dat_1)
  Seurat::DefaultAssay(pairedtag) <- "DNA"
  pairedtag <- Signac::RunTFIDF(pairedtag)
  pairedtag <-  Signac::FindTopFeatures(pairedtag, min.cutoff="q10")
  pairedtag <-  Signac::RunSVD(pairedtag)  # question: is this svd run with only the top variable features?
  set.seed(10)
  pairedtag <- Seurat::RunUMAP(pairedtag, 
                               reduction="lsi", 
                               dims=2:50,
                               reduction.name="umap.dna", 
                               reduction.key="dnaUMAP_")
  
  set.seed(10)
  pairedtag <- Seurat::FindMultiModalNeighbors(pairedtag, 
                                               reduction.list = list("pca", "lsi"), 
                                               dims.list = list(1:50, 2:50))
  set.seed(10)
  pairedtag <- Seurat::RunUMAP(pairedtag, 
                               nn.name = "weighted.nn", 
                               reduction.name = "wnn.umap", 
                               reduction.key = "wnnUMAP_")
  
  save(pairedtag, 
       file = paste0("../../../../data/Pairedtag_mousebrain_RNA-Histone/seurat_",
                     histone_names[i], ".RData"))
}