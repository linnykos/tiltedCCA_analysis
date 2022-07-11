# steady state plot
load("../../../out/main/10x_mouseembryo_steadystate.RData")
# load("../../out/main/10x_mouseembryo_steadystate.RData")

names(alignment_vec) <- colnames(mbrain)
alignment_vec1 <- alignment_vec[colnames(heatmap_mat1)]
df1 <- data.frame(y = alignment_vec1, x = 1:length(alignment_vec1))
# np_res1 <- npregfast::frfast(y ~ x, data = df1)
np_res1 <- stats::loess(y ~ x, data = df1)
col <- viridis::viridis(5)[3]

height <- 350
png(paste0("../../../out/figures/main/10x_mouseembryo_steadystate_forHeatmap1.png"),
    height = height, width = height*1800/360, units = "px", res = 500)
par(mar = c(0.1,1,0.1,0))
y <- np_res1$p[,1,1]; n <- length(y)
ylim <- quantile(alignment_vec1, probs = c(0.15,0.95))
ylim[1] <- min(ylim[1], min(y)); ylim[2] <- max(ylim[2], max(y))
tmp <- alignment_vec1; tol <- .01
tmp[which(tmp <= ylim[1]+tol)] <- NA
tmp[which(tmp >= ylim[2]-tol)] <- NA
plot(x = seq(0, 1, length = length(alignment_vec1)),
     y = tmp, 
     ylim = ylim,
     col = rgb(0.5, 0.5, 0.5, 0.1), main = "", 
     xaxt = "n", yaxt = "n", bty = "n", xlab = "", ylab = "",
     pch = 16)
lines(x = seq(0, 1, length = n), y = y, lwd = 9, col = "white")
lines(x = seq(0, 1, length = n), y = y, lwd = 6, col = col)
axis(2, labels = F, lwd = 2)
graphics.off()
