rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/Writeup14n/abseq_bm97Ref_tcca.RData")

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

bm[["distinctAB"]] <- Seurat::CreateAssayObject(data = t(multiSVD_obj$distinct_mat_2))

Seurat::DefaultAssay(bm) <- "distinctAB"
bm[["distinctAB"]]@var.features <- rownames(bm)
bm[["AB"]] <- NULL
bm[["RNA"]] <- NULL
bm[["umap"]] <- NULL
bm[["rna.umap"]] <- NULL
bm[["adt.umap"]] <- NULL
bm[["wnn.umap"]] <- NULL

tab <- table(bm$ct)
remove_celltypes <- names(tab)[which(tab <= 100)]
keep_idx <- rep(1, ncol(bm))
keep_idx[which(bm$ct %in% remove_celltypes)] <- 0
bm$keep <- keep_idx
bm <- subset(bm, keep == 1)
  
adt_distinct_de_list <- tiltedCCA:::differential_expression(seurat_obj = bm,
                                                            assay = "distinctAB",
                                                            idents = "ct",
                                                            test_use = "wilcox",
                                                            slot = "data")

save(adt_distinct_de_list, date_of_run, session_info,
     file = "../../../out/Writeup14n/abseq_bm97Ref_distinct_differential.RData")
