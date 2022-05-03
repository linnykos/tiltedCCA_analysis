rm(list=ls())
library(loomR)
library(Seurat)

# connect to the loom file
lfile <- loomR::connect(filename = "~/nzhanglab/data/10x_mouse_embryo/velocyto/10x_mouse_embryo.loom", 
                 mode = "r+", skip.validate = TRUE)
lfile
lfile[["col_attrs"]]
lfile[["col_attrs/CellID"]]
cellID <- lfile[["col_attrs/CellID"]][]
head(cellID)
cellID2 <- sapply(cellID, function(x){
  tmp <- strsplit(x, split = ":")[[1]][2]
  tmp <- substr(tmp, start = 0, stop = nchar(tmp)-1)
  paste0(tmp, "-1")
})
names(cellID2) <- NULL
head(cellID2)

# load in jane's annotations
mbrain.jane <- readRDS("~/project/tiltedCCA/data/10x_mouseembryo/data_tenx_labels_Jane.rds")
head(mbrain.jane@meta.data)

all(cellID2 %in% rownames(mbrain.jane@meta.data))
celltype_vec <- rep("NA", length(cellID2))
for(i in 1:length(cellID2)){
  cell_name <- cellID2[i]
  idx <- which(rownames(mbrain.jane@meta.data) == cell_name)
  celltype_vec[i] <- mbrain.jane$savercatLable[idx]
}
table(celltype_vec)

lfile$add.col.attribute(list(label.Savercat = celltype_vec), overwrite = TRUE)
lfile[["col_attrs"]]

x_vec <- lfile[["col_attrs/_X"]][]
y_vec <- lfile[["col_attrs/_Y"]][]

png("../../../out/figures/main/10x_mouseembryo_loompy-umap.png",
    height = 2500, width = 2500, units = "px", res = 300)
plot(x_vec, y_vec, pch = 16, col = as.numeric(as.factor(celltype_vec)))
graphics.off()

lfile$close_all()
