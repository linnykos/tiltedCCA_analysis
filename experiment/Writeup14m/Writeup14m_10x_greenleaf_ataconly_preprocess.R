rm(list=ls())
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)

# look at the gene activity matrix
Sys.setenv("VROOM_CONNECTION_SIZE" = 200000000)
mat <- readr::read_delim("~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_atac_counts.tsv.gz",
                         delim = "\t")
print(dim(mat))
save(mat,
     file = "~/nzhanglab/data/GSE162170_cortical_multiome/GSE162170_atac_counts.RData")
