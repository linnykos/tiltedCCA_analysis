rm(list=ls())
load("../../../out/main/10x_mouseembryo_preprocessed.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

n <- ncol(mbrain)
keep_vec <- rep(0, n)
keep_vec[which(mbrain$label_Savercat %in% c("Oligodendrocyte", "Radial glia", 
                                            "Forebrain GABAergic", "Neuroblast", 
                                            "Glioblast", "Cortical or hippocampal glutamatergic"))] <- 1
mbrain$keep <- keep_vec
mbrain <- subset(mbrain, keep == 1)

print("Working on RNA")
gene_de_list <- tiltedCCA:::differential_expression(seurat_obj = mbrain,
                                                    assay = "SCT",
                                                    idents = "label_Savercat",
                                                    test_use = "MAST",
                                                    slot = "counts")
save(gene_de_list, mbrain,
     date_of_run, session_info,
     file = "../../../out/main/10x_mouseembryo_differential.RData")

