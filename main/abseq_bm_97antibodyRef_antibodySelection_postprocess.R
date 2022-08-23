rm(list=ls())

library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../out/main/abseq_bm97Ref_varSelect.RData")

Seurat::DefaultAssay(bm) <- "AB"
for(ab_val in variable_selection_res$selected_variables){
  print(ab_val)
  
  plot1 <- Seurat::FeaturePlot(bm, 
                               features = ab_val,
                               reduction = "consensusUMAP",
                               slot = "scale.data")
  plot1 <- plot1 + Seurat::NoLegend() + Seurat::NoAxes()
  plot1 <- plot1 + ggplot2::ggtitle("")
  plot1 <- plot1 + ggplot2::theme(legend.text = ggplot2::element_text(size = 5))
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/abseq_bm97Ref_varSelect_consensusPCA-", ab_val, ".png"),
                  plot1, device = "png", width = 3, height = 3, units = "in",
                  dpi = 500)
}

mat <- bm[["AB"]]@scale.data

quantile(mat["CD324-AB",], probs = seq(0.9,1,length.out=11))
quantile(mat["CRACC-AB",], probs = seq(0.9,1,length.out=11))
