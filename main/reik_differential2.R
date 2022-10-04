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

df <- read.csv("~/project/tiltedCCA/data/mouse_cell_cycling/41467_2022_30545_MOESM5_ESM.txt", 
               header = F)

reik$celltype2 <- celltype
print("Working on RNA")
Seurat::DefaultAssay(reik) <- "RNA"
reik[["RNA"]]@var.features <- intersect(rownames(reik), df[,1])
gene_de_list <- tiltedCCA:::differential_expression(seurat_obj = reik,
                                                    assay = "RNA",
                                                    idents = "celltype2",
                                                    test_use = "MAST",
                                                    slot = "counts")
save(gene_de_list, 
     date_of_run, session_info,
     file = "../../../out/main/10x_reik_differential2.RData")

