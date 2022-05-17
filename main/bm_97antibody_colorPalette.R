col_palette <- c(
  "Late erythroid progenitor" = rgb(243, 159, 62, maxColorValue = 255),
  "Pre-B cells" = rgb(39, 108, 231, maxColorValue = 255),
  "Pro-B cells" = rgb(25, 72, 233, maxColorValue = 255),
  "Conventional dendritic cell 1" = rgb(167, 206, 163, maxColorValue = 255), #cDCs 1
  "Late promyelocytes" = rgb(109, 138, 40, maxColorValue = 255),
  "Megakaryocyte progenitors" = rgb(208, 50, 41, maxColorValue = 255),
  "Early promyelocytes" = rgb(188, 188, 120, maxColorValue = 255),
  "Lymphomyeloid prog" = rgb(74, 129, 61, maxColorValue = 255),
  "Early erythroid progenitor" = rgb(235, 110, 69, maxColorValue = 255),
  "Myelocytes" = rgb(89, 166, 69, maxColorValue = 255),
  "Conventional dendritic cell 2" = rgb(218, 226, 121, maxColorValue = 255), #cDCs 2
  "Eosinophil-basophil-mast cell progenitors" = rgb(161, 47, 46, maxColorValue = 255), #EoBaso progenitors 
  "Plasmacytoid dendritic cell progenitors" = rgb(25, 127, 29, maxColorValue = 255), 
  "HSCs & MPPs" = rgb(126, 26, 88, maxColorValue = 255), 
  "Plasma cells" = rgb(0, 0, 128, maxColorValue = 255), 
  "Erythro-myeloid progenitors" = rgb(205, 89, 136, maxColorValue = 255), 
  "CD4+ memory T cells" = rgb(164, 97, 186, maxColorValue = 255), 
  "Plasmacytoid dendritic cells" = rgb(31, 128, 28, maxColorValue = 255), 
  "Classical Monocytes" = rgb(54, 124, 113, maxColorValue = 255), 
  "NK cell progenitors" = rgb(181, 159, 196, maxColorValue = 255),  # only 11
  "Non-classical monocytes" = rgb(6, 199, 0, maxColorValue = 255),
  "CD8+CD103+ tissue resident memory T cells" = rgb(234, 117, 178, maxColorValue = 255), #CD8+CD103+ TRM T cells
  "CD56brightCD16- NK cells" = rgb(136, 107, 193, maxColorValue = 255), 
  "Class switched memory B cells" = rgb(68, 85, 130, maxColorValue = 255), 
  "CD56dimCD16+ NK cells" = rgb(57, 37, 142, maxColorValue = 255), 
  "CD8+ central memory T cells" = rgb(247, 198, 225, maxColorValue = 255), 
  "Pre-pro-B cells" = rgb(37, 79, 245, maxColorValue = 255),  # only 22
  "Nonswitched memory B cells" = rgb(38, 74, 153, maxColorValue = 255), 
  "CD8+ effector memory T cells" = rgb(231, 77, 153, maxColorValue = 255), 
  "CD8+ naive T cells" = rgb(228, 46, 118, maxColorValue = 255), 
  "Mature naive B cells" = rgb(150, 203, 252, maxColorValue = 255), 
  "GammaDelta T cells" = rgb(88, 52, 170, maxColorValue = 255), 
  "Immature B cells" = rgb(162, 233, 253, maxColorValue = 255),
  "CD4+ naive T cells" = rgb(94, 15, 150, maxColorValue = 255), 
  "CD69+PD-1+ memory CD4+ T cells" = rgb(234, 195, 240, maxColorValue = 255), # only 6
  "Small pre-B cell" = rgb(37, 37, 203, maxColorValue = 255), # only 12
  "Aberrant erythroid" = rgb(255, 205, 114, maxColorValue = 255), # only 2
  "NK T cells" = rgb(199, 188, 225, maxColorValue = 255) # only 1
)
col_palette <- col_palette[order(names(col_palette))]

# plot(1:length(col_palette), pch = 16, col = col_palette, cex = 5)

