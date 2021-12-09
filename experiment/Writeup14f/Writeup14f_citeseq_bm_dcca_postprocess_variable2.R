rm(list=ls())
load("../../../../out/Writeup14f/Writeup14f_citeseq_bm_dcca_postprocess_gene.RData")

for(i in 1:nrow(summary_mat2)){
  gene_name <- rownames(summary_mat2)[i]
  val <- quantile(sapply(1:length(adt_de_list), function(j){
    idx <- which(adt_de_list[[j]][,1] == gene_name)
    if(length(idx) == 0) return(1)
    adt_de_list[[j]][idx, "p_val"]
  }), probs = 0.5)
  
  summary_mat2[i,"p_val"] <- val
}

for(i in 1:nrow(summary_mat1)){
  if(i %% floor(nrow(summary_mat1)/10) == 0) cat('*')
  
  gene_name <- rownames(summary_mat1)[i]
  val <- quantile(sapply(1:length(gene_de_list), function(j){
    idx <- which(rownames(gene_de_list[[j]]) == gene_name)
    if(length(idx) == 0) return(1)
    gene_de_list[[j]][idx, "p_val"]
  }), probs = 0.5)
  
  summary_mat1[i,"p_val"] <- val
}

summary_mat1[,"p_val"] <- -log10(summary_mat1[,"p_val"])
summary_mat2[,"p_val"] <- -log10(summary_mat2[,"p_val"])

png("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_adt_exploration.png",
    height = 1250, width = 1250, res = 300, units = "px")
plot(NA,
     xlim = c(1, max(summary_mat2[,"p_val"])),
     ylim = c(0,1),
     bty = "n",
     xlab = "Separability (Average -log10(p-value))",
     ylab = "Alignment w/ common space (R^2)",
     main = "CITE-Seq: Bone marrow (Protein)",
     log = "x")
for(x in 10^(seq(0, 2.5, by = 0.5))){
  lines(rep(x,2), c(-1e5,1e5), lty = 2, col = "gray", lwd = 0.5)
}
for(y in seq(0, 1, by = 0.1)){
  lines(c(1,1e5), rep(y,2), lty = 2, col = "gray", lwd = 0.5)
}
points(summary_mat2[,3], summary_mat2[,1], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 1))
graphics.off()

########### ## ggrepel version of the plot
df <- as.data.frame(summary_mat2)
df$Name <- rownames(df)
p1 <- ggplot2::ggplot(df, ggplot2::aes(x = p_val, y = r_squared))
p1 <- p1 + ggplot2::geom_point()
p1 <- p1 + ggplot2::scale_x_log10()
p1 <- p1 + ggplot2::xlim(1, max(df$p_val))
p1 <- p1 + ggplot2::ylim(0, 1)
p1 <- p1 + ggrepel::geom_text_repel(data = df, ggplot2::aes(label = Name))
p1 <- p1 + Seurat::NoLegend()
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_adt_exploration2.png",
                p1, device = "png", width = 5, height = 5, units = "in")

###########################

summary_mat1[,"p_val"] <- pmax(summary_mat1[,"p_val"], 0.001)

png("../../../../out/figures/Writeup14f/Writeup14f_citeseq_bm_rna_exploration.png",
    height = 1250, width = 1250, res = 300, units = "px")
plot(NA,
     xlim = c(0.001, max(summary_mat1[,"p_val"])),
     ylim = c(0,1),
     bty = "n",
     xlab = "Separability (Average -log10(p-value))",
     ylab = "Alignment w/ common space (R^2)",
     main = "CITE-Seq: Bone marrow (Gene)",
     log = "x")
for(x in 10^(seq(-3, 2.5, by = 0.5))){
  lines(rep(x,2), c(-1e5,1e5), lty = 2, col = "gray", lwd = 0.5)
}
for(y in seq(0, 1, by = 0.1)){
  lines(c(0.0001,1e5), rep(y,2), lty = 2, col = "gray", lwd = 0.5)
}
points(summary_mat1[,3], summary_mat1[,1], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 0.2))
graphics.off()

