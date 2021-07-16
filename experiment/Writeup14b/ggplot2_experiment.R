rm(list=ls())

set.seed(10)
vec <- rnorm(200)
df <- data.frame(val = vec, idx = rank(vec), 
                 labeling = as.factor(c(rep(1,5),rep(0,195))),
                 text = paste0("p", 1:200))
df <- df[200:1,]

# see https://ggplot2.tidyverse.org/reference/scale_manual.html
# see https://stackoverflow.com/questions/15706281/controlling-order-of-points-in-ggplot2-in-r
# see https://stackoverflow.com/questions/34588232/subset-parameter-in-layers-is-no-longer-working-with-ggplot2-2-0-0
# see https://cran.r-project.org/web/packages/ggrepel/vignettes/ggrepel.html
# see https://stackoverflow.com/questions/15624656/label-points-in-geom-point
# see https://mran.microsoft.com/snapshot/2017-08-20/web/packages/ggrepel/vignettes/ggrepel.html
# see https://ggrepel.slowkow.com/articles/examples.html for defaults
p1 <- ggplot2::ggplot(df, ggplot2::aes(x = idx, y = val)) 
p1 <- p1 + ggplot2::geom_point(ggplot2::aes(color = labeling))
# p1 <- p1 + ggplot2::geom_point(data = subset(df, labeling == 1), ggplot2::aes(color = labeling))
p1 <- p1 + ggplot2::scale_colour_manual(values=c("black", "red"))
p1 <- p1 + ggrepel::geom_text_repel(data = subset(df, labeling == 1), ggplot2::aes(label = text, color = labeling),
                                    box.padding = ggplot2::unit(0.5, 'lines'),
                                    point.padding = ggplot2::unit(1.6, 'lines'),
                                    nudge_y = 1)
p1 <- p1 + Seurat::NoLegend()
p1


