rm(list=ls())
load("../../../out/main/citeseq_pbmc224_preprocessed.RData")

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

pbmc$celltype.l2_custom <- celltype_l3
print("Working on RNA")
gene_de_list <- tiltedCCA:::differential_expression(seurat_obj = pbmc,
                                                    assay = "SCT",
                                                    idents = "celltype.l2_custom",
                                                    test_use = "MAST",
                                                    slot = "counts")
save(gene_de_list, pbmc,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_differential.RData")

########

pbmc[["ADT"]]@var.features <- rownames(pbmc[["ADT"]])
print("Working on ADT")
adt_de_list <- tiltedCCA:::differential_expression(seurat_obj = pbmc,
                                                    assay = "ADT",
                                                    idents = "celltype.l2_custom",
                                                    test_use = "wilcox",
                                                    slot = "data")

save(gene_de_list, adt_de_list, pbmc,
     date_of_run, session_info,
     file = "../../../out/main/citeseq_pbmc224_differential.RData")
