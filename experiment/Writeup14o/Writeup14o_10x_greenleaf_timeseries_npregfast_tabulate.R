rm(list=ls())
library(Seurat)
library(Signac)
library(tiltedCCA)

load("../../../../out/Writeup14o/Writeup14o_10x_greenleaf_timeseries_npregfast.RData")

p <- ncol(rna_common)
sd_vec <- sapply(1:p, function(j){
  
})