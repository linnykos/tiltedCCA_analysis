rm(list=ls())
load("../../../out/main/citeseq_pbmc224_varSelect_enrichment.RData")
source("pbmc_224antibody_colorPalette.R")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

png("../../../out/figures/main/citeseq_pbmc224_varSelect_enrichment.png",
    height = 2000, width = 3500, units = "px", res = 300)
par(mfrow = c(1,4), mar = c(4,4,4,0.5))
barplot(enrichment_selected$enrichment$df[,"value"],
        col = col_palette[enrichment_selected$enrichment$df[,"celltype"]],
        ylim = c(0,1),
        names.arg = enrichment_selected$enrichment$df[,"celltype"],
        cex.names = 0.5,
        main = "Selected")

barplot(enrichment_alt_1$enrichment$df[,"value"],
        col = col_palette[enrichment_alt_1$enrichment$df[,"celltype"]],
        ylim = c(0,1),
        names.arg = enrichment_alt_1$enrichment$df[,"celltype"],
        cex.names = 0.5,
        main = "Highest p-value")

barplot(enrichment_alt_2$enrichment$df[,"value"],
        col = col_palette[enrichment_alt_2$enrichment$df[,"celltype"]],
        ylim = c(0,1),
        names.arg = enrichment_alt_2$enrichment$df[,"celltype"],
        cex.names = 0.5,
        main = "Least R-squared")

barplot(enrichment_alt_3$enrichment$df[,"value"],
        col = col_palette[enrichment_alt_3$enrichment$df[,"celltype"]],
        ylim = c(0,1),
        names.arg = enrichment_alt_3$enrichment$df[,"celltype"],
        cex.names = 0.5,
        main = "All RNA")

graphics.off()

mean(enrichment_selected$enrichment$df[,"value"])
mean(enrichment_alt_1$enrichment$df[,"value"])
mean(enrichment_alt_2$enrichment$df[,"value"])
mean(enrichment_alt_3$enrichment$df[,"value"])

mean(enrichment_selected$enrichment$df[,"sd"])
mean(enrichment_alt_1$enrichment$df[,"sd"])
mean(enrichment_alt_2$enrichment$df[,"sd"])
mean(enrichment_alt_3$enrichment$df[,"sd"])

median(enrichment_selected$enrichment$df[,"value"])
median(enrichment_alt_1$enrichment$df[,"value"])
median(enrichment_alt_2$enrichment$df[,"value"])
median(enrichment_alt_3$enrichment$df[,"value"])