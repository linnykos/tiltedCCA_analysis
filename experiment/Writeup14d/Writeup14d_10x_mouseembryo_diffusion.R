rm(list=ls())
load("../../../../out/Writeup14d/Writeup14d_10x_mouseembryo_dcca.RData")

cell_idx <- which(as.character(mbrain2@meta.data$wsnn_res.2) %in% c("17", "16", "0", "4", "2", "1", "5", "7", "12", "20"))
cell_idx <- cell_idx[-which(!mbrain2@meta.data$label_Savercat[cell_idx] %in% c("Radial glia", "Neuroblast", "Cortical or hippocampal glutamatergic"))]

obj <- 