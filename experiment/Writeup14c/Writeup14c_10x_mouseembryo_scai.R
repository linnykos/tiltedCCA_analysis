rm(list=ls()); set.seed(10)

library(Seurat); library(Signac); library(scAI)
library(dplyr); library(cowplot)

load("../../../../out/Writeup14/Writeup14_10x_mouseembryo_preprocess.RData")

metadata <- mbrain@meta.data
cell_idx <- which(metadata$label_Savercat %in% c("Oligodendrocyte", "Radial glia", 
                                                 "Forebrain GABAergic", "Neuroblast", 
                                                 "Glioblast", "Cortical or hippocampal glutamatergic"))

X <- list(mbrain[["SCT"]]@data[,cell_idx], 
          mbrain[["ATAC"]]@data[Seurat::VariableFeatures(object = mbrain),cell_idx])
df <- data.frame(celltype = mbrain@meta.data$label_Savercat[cell_idx])

rm(list = c("mbrain"))

scAI_outs <- scAI::create_scAIobject(raw.data = X, do.sparse = T)
rm(list = c("X"))
scAI_outs <- scAI::preprocessing(scAI_outs, assay = NULL)
scAI_outs <- scAI::addpData(scAI_outs, pdata = df, pdata.name = "Cell types")

set.seed(10)
scAI_outs <- scAI::run_scAI(scAI_outs, K = 30, nrun = 1, do.fast = T)

set.seed(10)
scAI_outs <- scAI::reducedDims(scAI_outs, method = "umap")

save.image(file = "../../../../out/Writeup14c/10x_mouseembryo_scai.RData")

