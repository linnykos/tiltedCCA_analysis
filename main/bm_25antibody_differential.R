rm(list=ls())
load("../../../out/main/citeseq_bm25_preprocessed.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

celltype_l3 <- as.character(bm$celltype.l2)
celltype_l3[celltype_l3 %in% c("NK", "CD56 bright NK")] <- "NK_all"
celltype_l3[celltype_l3 %in% c("MAIT", "gdT")] <- "MAIT-gdT"
celltype_l3[celltype_l3 %in% c("Memory B", "Naive B")] <- "B"
celltype_l3[celltype_l3 %in% c("Prog_B 1", "Prog_B 2")] <- "Prog_B"
celltype_l3[celltype_l3 %in% c("Prog_Mk", "Plasmablast", "LMPP", "Treg")] <- NA

bm$celltype.l3 <- celltype_l3
gene_de_list <- tiltedCCA:::differential_expression(seurat_obj = bm,
                                                    assay = "RNA",
                                                    idents = "celltype.l3",
                                                    test_use = "MAST",
                                                    slot = "counts")
save(gene_de_list, bm,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_bm25_differential.RData")

########

adt_de_list <- tiltedCCA:::differential_expression(seurat_obj = bm,
                                                    assay = "ADT",
                                                    idents = "celltype.l3",
                                                    test_use = "wilcox",
                                                    slot = "data")

save(gene_de_list, adt_de_list, bm,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_bm25_differential.RData")
