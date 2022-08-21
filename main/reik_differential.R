rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/10x_reik_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

celltype <- as.character(reik$celltype)
rm_celltypes <- names(table(celltype))[which(table(celltype) <= 200)]
celltype[celltype %in% rm_celltypes] <- NA

reik$celltype2 <- celltype
print("Working on RNA")
Seurat::DefaultAssay(reik) <- "RNA"
gene_de_list <- tiltedCCA:::differential_expression(seurat_obj = reik,
                                                    assay = "RNA",
                                                    idents = "celltype2",
                                                    test_use = "MAST",
                                                    slot = "counts")
save(gene_de_list, reik,
     date_of_run, session_info,
     file = "../../../out/main/10x_reik_differential.RData")

