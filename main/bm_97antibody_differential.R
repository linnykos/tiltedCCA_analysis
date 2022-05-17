rm(list=ls())
load("../../../out/main/abseq_bm97_preprocessed.RData")

library(Seurat)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(1, ncol(bm))
keep_vec[which(bm$ct %in% c("Aberrant erythroid", "CD69+PD-1+ memory CD4+ T cells",
                      "NK cell progenitors", "NK T cells",
                      "Plasma cells", "Pre-pro-B cells", "Small pre-B cell"))] <- 0
bm$keep <- keep_vec
bm <- subset(bm, keep == 1)
bm$ct <- droplevels(bm$ct)

ct2 <- as.character(bm$ct)
ct2[ct2 %in% c("Conventional dendritic cell 1", "Conventional dendritic cell 2")] <- "cDC"

bm$ct2 <- ct2
print("Working on RNA")
gene_de_list <- tiltedCCA:::differential_expression(seurat_obj = bm,
                                                    assay = "RNA",
                                                    idents = "ct2",
                                                    test_use = "MAST",
                                                    slot = "counts")
save(gene_de_list, bm,
     date_of_run, session_info,
     file = "../../../out/main/abseq_bm97_differential.RData")

########

bm[["AB"]]@var.features <- rownames(bm[["AB"]])
print("Working on ADT")
adt_de_list <- tiltedCCA:::differential_expression(seurat_obj = bm,
                                                   assay = "AB",
                                                   idents = "ct2",
                                                   test_use = "wilcox",
                                                   slot = "data")

save(gene_de_list, adt_de_list, bm,
     date_of_run, session_info,
     file = "../../../out/main/abseq_bm97_differential.RData")
