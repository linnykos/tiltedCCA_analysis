plot_custom_cluster <- function(dimred,
                                membership_vec,
                                color_vec,
                                filename,
                                main,
                                height = 1200,
                                width = 1200,
                                xlab = "PC1",
                                ylab = "PC2",
                                ...){
  png(filename, height = height, width = width, units = "px", res = 300)
  graphics::plot(dimred[,1], dimred[,2],
                 asp = T,
                 pch = 16,
                 col = color_vec[membership_vec],
                 main = main,
                 xlab = xlab, ylab = ylab, ...)
  graphics.off()
}