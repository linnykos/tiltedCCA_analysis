rm(list=ls())
library(Seurat)
library(tiltedCCA)

load("../../../out/main/abseq_bm97Ref_preprocessed.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

keep_vec <- rep(1, ncol(bm))
tab <- table(bm$ct)
keep_vec[which(bm$ct %in% names(tab)[which(tab <= 150)])] <- 0
bm$keep <- keep_vec
bm <- subset(bm, keep == 1)
bm$ct <- droplevels(factor(bm$ct))
print(paste0("Number of levels: ", length(levels(bm$ct))))

print("Working on RNA")
gene_de_list <- tiltedCCA:::differential_expression(seurat_obj = bm,
                                                    assay = "RNA",
                                                    idents = "ct",
                                                    test_use = "MAST",
                                                    slot = "counts")
save(gene_de_list, bm,
     date_of_run, session_info,
     file = "../../../out/main/abseq_bm97Ref_differential.RData")

########

bm[["AB"]]@var.features <- rownames(bm[["AB"]])
print("Working on ADT")
adt_de_list <- tiltedCCA:::differential_expression(seurat_obj = bm,
                                                   assay = "AB",
                                                   idents = "ct",
                                                   test_use = "wilcox",
                                                   slot = "data")

save(gene_de_list, adt_de_list, bm,
     date_of_run, session_info,
     file = "../../../out/main/abseq_bm97Ref_differential.RData")
