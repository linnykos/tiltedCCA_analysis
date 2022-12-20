col_palette <- c(
  "Late erythroid progenitor" = rgb(243, 159, 62, maxColorValue = 255),
  "Pre-B cells" = rgb(39, 108, 231, maxColorValue = 255),
  "Pro-B cells" = rgb(39, 108, 231, maxColorValue = 255),
  "Conventional dendritic cell 1" = rgb(89, 166, 69, maxColorValue = 255), #cDCs 1
  "Late promyelocytes" = rgb(89, 166, 69, maxColorValue = 255),
  "Megakaryocyte progenitors" = rgb(243, 159, 62, maxColorValue = 255),
  "Early promyelocytes" = rgb(89, 166, 69, maxColorValue = 255),
  "Lymphomyeloid prog" = rgb(89, 166, 69, maxColorValue = 255),
  "Early erythroid progenitor" = rgb(243, 159, 62, maxColorValue = 255),
  "Myelocytes" = rgb(89, 166, 69, maxColorValue = 255),
  "Conventional dendritic cell 2" = rgb(89, 166, 69, maxColorValue = 255), #cDCs 2
  "Eosinophil-basophil-mast cell progenitors" = rgb(243, 159, 62, maxColorValue = 255), #EoBaso progenitors 
  "Plasmacytoid dendritic cell progenitors" = rgb(89, 166, 69, maxColorValue = 255), 
  "HSCs & MPPs" = rgb(243, 159, 62, maxColorValue = 255), 
  "Plasma cells" = rgb(164, 97, 186, maxColorValue = 255), 
  "Erythro-myeloid progenitors" = rgb(243, 159, 62, maxColorValue = 255),
  "CD4+ cytotoxic T cells" = rgb(164, 97, 186, maxColorValue = 255), 
  "CD4+ memory T cells" = rgb(164, 97, 186, maxColorValue = 255), 
  "Plasmacytoid dendritic cells" = rgb(89, 166, 69, maxColorValue = 255),
  "Classical Monocytes" = rgb(89, 166, 69, maxColorValue = 255),
  "NK cell progenitors" = rgb(181, 159, 196, maxColorValue = 255),  # only 11
  "Non-classical monocytes" = rgb(89, 166, 69, maxColorValue = 255),
  "CD8+CD103+ tissue resident memory T cells" = rgb(234, 117, 178, maxColorValue = 255), #CD8+CD103+ TRM T cells
  "CD56brightCD16- NK cells" = rgb(181, 159, 196, maxColorValue = 255),
  "Class switched memory B cells" = rgb(39, 108, 231, maxColorValue = 255),
  "CD56dimCD16+ NK cells" = rgb(181, 159, 196, maxColorValue = 255),
  "CD8+ central memory T cells" = rgb(234, 117, 178, maxColorValue = 255),
  "Pre-pro-B cells" = rgb(39, 108, 231, maxColorValue = 255),  # only 22
  "Nonswitched memory B cells" = rgb(39, 108, 231, maxColorValue = 255),
  "CD8+ effector memory T cells" = rgb(234, 117, 178, maxColorValue = 255),
  "CD8+ naive T cells" = rgb(234, 117, 178, maxColorValue = 255),
  "Mature naive B cells" = rgb(39, 108, 231, maxColorValue = 255),
  "GammaDelta T cells" = rgb(164, 97, 186, maxColorValue = 255), 
  "Immature B cells" = rgb(39, 108, 231, maxColorValue = 255),
  "CD4+ naive T cells" = rgb(164, 97, 186, maxColorValue = 255), 
  "CD69+PD-1+ memory CD4+ T cells" = rgb(164, 97, 186, maxColorValue = 255),  # only 6
  "Small pre-B cell" = rgb(39, 108, 231, maxColorValue = 255), # only 12
  "Aberrant erythroid" = rgb(243, 159, 62, maxColorValue = 255), # only 2
  "NK T cells" = rgb(181, 159, 196, maxColorValue = 255), # only 1
  "CD11c+ memory B cells Monocyte-like blasts" = rgb(39, 108, 231, maxColorValue = 255),
  "Mesenchymal cells_1" = rgb(89, 166, 69, maxColorValue = 255),
  "Mesenchymal cells_2" = rgb(89, 166, 69, maxColorValue = 255)
)
col_palette <- col_palette[order(names(col_palette))]

plot(1:length(col_palette), pch = 16, col = col_palette, cex = 5)

