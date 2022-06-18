rm(list=ls())
load("../../../out/main/citeseq_pbmc224_tiltedcca.RData")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

celltype_l3 <- as.character(pbmc$celltype.l2)
celltype_l3[celltype_l3 %in% c("NK", "NK Proliferating", "NK_CD56bright")] <- "NK_all"
celltype_l3[celltype_l3 %in% c("B intermediate", "B memory")] <- "B_memory_intermediate"
celltype_l3[celltype_l3 %in% c("ILC", "ASDC", "CD4 Proliferating", "CD8 Proliferating", "cDC1", "dnT", "Eryth", "HSPC", "Plasmablast")] <- NA
pbmc$celltype.l3 <- celltype_l3

pbmc[["distinctADT"]] <- Seurat::CreateAssayObject(data = t(multiSVD_obj$distinct_mat_2))

print("Working on RNA")
adt_distinct_de_list <- tiltedCCA:::differential_expression(seurat_obj = pbmc,
                                                    assay = "distinctADT",
                                                    idents = "celltype.l3",
                                                    test_use = "wilcox",
                                                    slot = "data")
save(adt_distinct_de_list, pbmc,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_differential_distinct.RData")