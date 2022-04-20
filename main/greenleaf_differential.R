rm(list=ls())
load("../../../out/main/10x_greenleaf_preprocessed.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

celltype_l2 <- as.character(greenleaf$celltype)
celltype_l2[celltype_l2 %in% c("EC/Peric.", "SP")] <- NA

greenleaf$celltype.l2 <- celltype_l2
print("Working on RNA")
gene_de_list <- tiltedCCA:::differential_expression(seurat_obj = greenleaf,
                                                    assay = "SCT",
                                                    idents = "celltype.l2",
                                                    test_use = "MAST",
                                                    slot = "counts")
save(gene_de_list, greenleaf,
     date_of_run, session_info,
     file = "../../../out/main/10x_greenleaf_differential.RData")

