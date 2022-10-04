rm(list=ls())
load("../../../out/main/citeseq_bm25_tcca.RData")
source("bm_25antibody_colorPalette.R")

library(Seurat)
library(Signac)
library(tiltedCCA)

set.seed(10)
date_of_run <- Sys.time()
session_info <- devtools::session_info()

membership_vec <- factor(bm$celltype.l2)

cell_enrichment_res <- tiltedCCA:::postprocess_cell_enrichment(
  input_obj = multiSVD_obj, 
  membership_vec = membership_vec, 
  max_subsample = 1000,
  verbose = 1
)

cbind(cell_enrichment_res$enrichment_common$df[,1], round(cell_enrichment_res$enrichment_common$df[,2],2))
cbind(cell_enrichment_res$enrichment_distinct_1$df[,1], round(cell_enrichment_res$enrichment_distinct_1$df[,2],2))
cbind(cell_enrichment_res$enrichment_distinct_2$df[,1], round(cell_enrichment_res$enrichment_distinct_2$df[,2],2))

cbind(cell_enrichment_res$enrichment_common$df[,c(1:2)], cell_enrichment_res$enrichment_distinct_1$df[,2])
cbind(cell_enrichment_res$enrichment_common$df[,c(1:2)], cell_enrichment_res$enrichment_distinct_2$df[,2])

save(cell_enrichment_res, 
     date_of_run, session_info,
     file = "../../../out/main/citeseq_bm25_cell_enrichment.RData")

##########################

# make a bar plot
source("../tiltedCCA_analysis/main/bm_25antibody_colorPalette.R")
load("../../out/main/citeseq_bm25_cell_enrichment.RData")

celltype_list <- list(Progenitors = c("GMP", "HSC", "LMPP", "Prog_DC", "Prog_Mk", "Prog_RBC", "Prog_B 1", "Prog_B 2"),
                      CD4 = c("CD4 Memory", "CD4 Naive", "Treg"),
                      CD8 = c("CD8 Effector_1", "CD8 Effector_2", "CD8 Memory_1", "CD8 Memory_2", "CD8 Naive"))
celltype_all <- unlist(celltype_list)
max_val <- max(max(cell_enrichment_res$enrichment_distinct_1$df[which(cell_enrichment_res$enrichment_distinct_1$df[,1] %in% celltype_all), 2]),
               max(cell_enrichment_res$enrichment_distinct_2$df[which(cell_enrichment_res$enrichment_distinct_2$df[,1] %in% celltype_all), 2]))
max_val <- ceiling(max_val * 10)/10

for(k in 1:2){
  obj_name <- ifelse(k == 1, "enrichment_distinct_1", "enrichment_distinct_2")
  
  vec <- unlist(lapply(1:3, function(i){
    idx_vec <- sapply(celltype_list[[i]], function(x){which(cell_enrichment_res[[obj_name]]$df[,1] == x)})
    tmp <- cell_enrichment_res[[obj_name]]$df[idx_vec, 2]
    names(tmp) <- celltype_list[[i]]
    tmp
  }))
  
  png(paste0("../../out/figures/main/citeseq_bm25_cell_enrichment_distinct-", k, "_cleaned.png"),
      height = 3500, width = 2000, units = "px", res = 500)
  par(mar = c(1.2,4,0,0))
  barplot(vec,  ylim = c(0, max_val), 
          space = c(rep(0,length(celltype_list[[1]])),
                   2, rep(0,length(celltype_list[[2]])-1),
                   2, rep(0,length(celltype_list[[3]])-1)),
          col = col_palette[names(vec)], 
          names.arg = rep("", length(vec)),
          xaxt = "n", yaxt = "n", bty = "n")
  graphics::axis(2,
                 cex.axis = 2,
                 lwd = 2,
                 lwd.ticks = 2)
  graphics.off()
}

######################

y_vec <- cell_enrichment_res$enrichment_common$df[,"value"]
names(y_vec) <- cell_enrichment_res$enrichment_common$df[,"celltype"]
x_vec <- cell_enrichment_res$enrichment_distinct_2$df[,"value"]
names(x_vec) <- cell_enrichment_res$enrichment_distinct_2$df[,"celltype"]

png("../../../out/figures/main/bm_25antibody_cell_enrichment_common-vs-distinct2.png",
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(4, 4, 0.5, 0.5))
plot(NA, xlim = c(0,1), ylim = c(0,1),
     main = "", xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n", asp = T)
for(x in seq(0,1,by=0.1)){
  lines(rep(x,2), c(-10,10), col = "gray", lty = 3, lwd = 2)
}
for(y in seq(0,1,by=0.1)){
  lines(c(-10,10), rep(y,2), col = "gray", lty = 3, lwd = 2)
}
lines(c(-10,10), c(-10,10), col = 2, lty = 2, lwd = 3)

col_vec <- col_palette[names(x_vec)]
points(x_vec, y_vec, col = 1, cex = 4, pch = 16)
points(x_vec, y_vec, col = "white", cex = 3, pch = 16)
points(x_vec, y_vec, col = col_vec, cex = 2.5, pch = 16)

axis(1, cex.axis = 1.25, cex.lab = 1.25)
axis(2, cex.axis = 1.25, cex.lab = 1.25)

graphics.off()


y_vec <- cell_enrichment_res$enrichment_common$df[,"value"]
names(y_vec) <- cell_enrichment_res$enrichment_common$df[,"celltype"]
x_vec <- -cell_enrichment_res$enrichment_distinct_1$df[,"value"]
names(x_vec) <- cell_enrichment_res$enrichment_distinct_1$df[,"celltype"]

png("../../../out/figures/main/bm_25antibody_cell_enrichment_common-vs-distinct1.png",
    height = 2500, width = 2500, res = 500, units = "px")
par(mar = c(4, 4, 0.5, 0.5))
plot(NA, xlim = c(-1,0), ylim = c(0,1),
     main = "", xlab = "", ylab = "",
     xaxt = "n", yaxt = "n", bty = "n", asp = T)
for(x in seq(-1,0,by=0.1)){
  lines(rep(x,2), c(-10,10), col = "gray", lty = 3, lwd = 2)
}
for(y in seq(0,1,by=0.1)){
  lines(c(-10,10), rep(y,2), col = "gray", lty = 3, lwd = 2)
}
lines(c(-10,10), c(10,-10), col = 2, lty = 2, lwd = 3)

col_vec <- col_palette[names(x_vec)]
points(x_vec, y_vec, col = 1, cex = 4, pch = 16)
points(x_vec, y_vec, col = "white", cex = 3, pch = 16)
points(x_vec, y_vec, col = col_vec, cex = 2.5, pch = 16)

axis(1, cex.axis = 1.25, cex.lab = 1.25)
# axis(2, cex.axis = 1.25, cex.lab = 1.25)
graphics.off()

