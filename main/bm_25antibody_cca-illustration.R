rm(list=ls())
load("../../../out/main/citeseq_bm25_tcca.RData")
source("bm_25antibody_colorPalette.R")

celltype_vec <- bm$celltype.l2
col_vec <- sapply(celltype_vec, function(celltype){
  col_palette[which(names(col_palette) == celltype)]
})

vec1 <- multiSVD_obj$cca_obj$score_1[,1]
vec2 <- multiSVD_obj$cca_obj$score_2[,1]

png("../../../out/figures/main/citeseq_bm25_cca_leading-scores.png",
    height = 1200, width = 1200, units = "px", res = 500)
par(mar = c(0.5, 4, 0.5, 0.5), bg = NA)
plot(NA, ylim = range(c(vec1, vec2)), xlim = c(-.5, 1.5), 
     ylab = "Leading canonical score", xlab = "", bty = "n", xaxt = "n",
     yaxt = "n")
axis(side = 2, labels = F)
points(y = vec1, x = rep(0, length(vec1)) + runif(length(vec1), min = -.3, max = .3),
       pch = 16, col = col_vec, cex = 0.25)
points(y = vec2, x = rep(1, length(vec2)) + runif(length(vec1), min = -.3, max = .3), 
       pch = 16, col = col_vec, cex = 0.25)
graphics.off()