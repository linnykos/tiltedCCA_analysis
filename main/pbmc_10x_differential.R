rm(list=ls())
load("../../../out/main/10x_pbmc_preprocessed.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

celltype <- as.character(pbmc$predicted.id)
celltype[celltype %in% c("NK", "NK Proliferating", "NK_CD56bright")] <- "NK_all"
celltype[celltype %in% c("B naive", "B intermediate", "B memory")] <- "B"
celltype[celltype %in% c("CD8 TCM", "CD8 TEM")] <- "CD8_memory"
celltype[celltype %in% c("ILC", "Plasmablast", "HSPC", "pDC", "ASDC", "CD4 CTL", 
                         "CD4 Proliferating", "CD4 TEM", "cDC1", "cDC2", "dnT",
                         "gdT", "ILC", "Platelet")] <- NA

pbmc$celltype <- celltype
print("Working on RNA")
gene_de_list <- tiltedCCA:::differential_expression(seurat_obj = pbmc,
                                                    assay = "SCT",
                                                    idents = "celltype",
                                                    test_use = "MAST",
                                                    slot = "counts")
save(gene_de_list, pbmc,
     date_of_run, session_info,
     file = "../../../out/main/10x_pbmc_differential.RData")

