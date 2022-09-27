rm(list=ls())
library(Seurat)
library(Signac)

load("../../../out/main/10x_greenleaf_preprocessed.RData")
source("greenleaf_colorPalette.R")

keep_vec <- rep(1, ncol(greenleaf))
keep_vec[greenleaf$celltype %in% c("EC/Peric.", "IN1", "IN2", "IN3", "mGPC/OPC", "SP")] <- 0
greenleaf$keep <- keep_vec
greenleaf <- subset(greenleaf, keep == 1)

Seurat::Idents(greenleaf) <- "celltype"

genes <- c("NTRK2", "CNTN1", "GRM1", "SEMA3E",
           "ASCL1", "EOMES", "PAX6", "ZNF521",
           "NEUROD1", "KCNH8", "MEIS2", "WNT11")

for(gene in genes){
  plot1 <- Signac::CoveragePlot(
    object = greenleaf,
    region = gene,
    extend.upstream = 5000,
    extend.downstream = 5000
  )
  ggplot2::ggsave(filename = paste0("../../../out/figures/main/10x_greenleaf_coverage_", gene, ".png"),
                  plot1, device = "png", width = 6, height = 5, units = "in")
}
