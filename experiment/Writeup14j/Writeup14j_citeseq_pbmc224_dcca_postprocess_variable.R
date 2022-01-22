rm(list=ls())
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca.RData")
load("../../../../out/Writeup14i/Writeup14i_citeseq_pbmc224_dcca_variable_full_de.RData")
source("../Writeup14f/gene_exploration.R")
load("../../../../out/Writeup14j/Writeup14j_citeseq_pbmc224_dcca_tmp.RData")
dcca_res <- res
class(dcca_res) <- "dcca"

dcca_decomp <- multiomicCCA::dcca_decomposition(dcca_res)

summary_mat <- compute_variable_summary(mat = mat_2b, 
                                        common_mat = dcca_decomp$common_mat_2,
                                        metacell_clustering = NA)
for(i in 1:nrow(summary_mat)){
  gene_name <- rownames(summary_mat)[i]
  
  val <- stats::median(sapply(1:length(de_list), function(j){
    idx <- which(rownames(de_list[[j]]) == gene_name)
    if(length(idx) == 0) return(1)
    de_list[[j]][idx, "p_val"]
  }))
  
  summary_mat[i,"kl_div"] <- val
}
summary_mat[,"kl_div"] <- -log10(summary_mat[,"kl_div"])

png("../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_adt_exploration.png",
    height = 1250, width = 1250, res = 300, units = "px")
tmp <- summary_mat[,"kl_div"]; tmp <- tmp[tmp > 0]
xmin <- floor(log10(min(tmp)))
summary_mat[which(summary_mat[,"kl_div"] == 0),"kl_div"] <- 10^xmin
xmax <- ceiling(max(log10(summary_mat[,"kl_div"])))
plot(NA,
     xlim = c(xmin, xmax),
     ylim = c(0,1),
     bty = "n",
     xaxt = "n",
     xlab = "Separability (Median -Log10(p-value))",
     ylab = "Alignment w/ common space (R^2)",
     main = "CITE-Seq: PBMC (Protein)")
tmp <- seq(xmin, xmax, by = 1)
xaxt_vec <- sort(unique(unlist(sapply(1:(length(tmp)-1), function(i){
  lower <- 10^(tmp[i]); upper <- 10^(tmp[i+1])
  seq_vec <- seq(lower, upper, length.out = 10)
  log10(seq_vec)
}))))
axis(1, at = xaxt_vec, labels = rep("", length(xaxt_vec)))
axis(1, at = tmp, labels = as.character(10^tmp))
for(y_val in seq(0,1,by=0.1)){
  lines(c(-1e5,1e5), rep(y_val,2), lty = 2, col = "gray", lwd = 1)
}
tmp <- pmax(seq(0,length(xaxt_vec), by = 2),1)
major <- c(1, tmp[tmp %% 10 == 0])
minor <- tmp[!tmp %in% major]
for(x_val in xaxt_vec[minor]){
  lines(rep(x_val,2), c(-1e5,1e5), lty = 3, col = "gray", lwd = 0.5)
}
for(x_val in xaxt_vec[major]){
  lines(rep(x_val,2), c(-1e5,1e5), lty = 2, col = "gray", lwd = 1.5)
}

points(log10(summary_mat[,"kl_div"]), summary_mat[,"r_squared"], pch = 16,
       col = rgb(0.5, 0.5, 0.5, 1))
graphics.off()


########### ## ggrepel version of the plot
df <- as.data.frame(summary_mat)
df$Name <- rownames(df)
p1 <- ggplot2::ggplot(df, ggplot2::aes(x = kl_div, y = r_squared))
p1 <- p1 + ggplot2::geom_point()
p1 <- p1 + ggplot2::scale_x_log10()
# p1 <- p1 + ggplot2::xlim(1, max(df$kl_div))
p1 <- p1 + ggplot2::ylim(0, 1)
p1 <- p1 + ggrepel::geom_text_repel(data = df, ggplot2::aes(label = Name))
p1 <- p1 + Seurat::NoLegend()
ggplot2::ggsave(filename = "../../../../out/figures/Writeup14j/Writeup14j_citeseq_pbmc224_adt_exploration2.png",
                p1, device = "png", width = 5, height = 5, units = "in")

