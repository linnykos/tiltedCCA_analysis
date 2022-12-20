col_palette <- c(
  "CD14 Mono" = rgb(255, 108, 145, maxColorValue = 255),
  "CD16 Mono" = rgb(255, 104, 159, maxColorValue = 255),
  "CD4 Memory" = rgb(231, 134, 26, maxColorValue = 255),
  "CD4 Naive"= rgb(224, 139, 0, maxColorValue = 255),
  "CD56 bright NK" = rgb(0.7, 0.7, 0.7),
  "CD8 Effector_1" = rgb(2, 190, 108, maxColorValue = 255),
  "CD8 Effector_2" = rgb(0, 188, 89, maxColorValue = 255),
  "CD8 Memory_1" = rgb(57, 182, 0, maxColorValue = 255),
  "CD8 Memory_2" = rgb(0, 186, 66, maxColorValue = 255),
  "CD8 Naive" = rgb(1, 184, 31, maxColorValue = 255),
  "cDC2" = rgb(255, 98, 188, maxColorValue = 255),
  "gdT" = rgb(133, 173, 0, maxColorValue = 255),
  "GMP" = rgb(255, 97, 201, maxColorValue = 255),
  "HSC" = rgb(247, 99, 224, maxColorValue = 255),
  "LMPP" = rgb(240, 102, 234, maxColorValue = 255),
  "MAIT" = rgb(149, 169, 0, maxColorValue = 255),
  "Memory B" = rgb(0, 180, 239, maxColorValue = 255),
  "Naive B" = rgb(2, 165, 255, maxColorValue = 255),
  "NK" = rgb(0.6, 0.6, 0.6),
  "pDC" = rgb(255, 101, 174, maxColorValue = 255),
  "Plasmablast" = rgb(1, 191, 196, maxColorValue = 255),
  "Prog_B 1" = rgb(172, 136, 255, maxColorValue = 255),
  "Prog_B 2" = rgb(121, 151, 255, maxColorValue = 255),
  "Prog_DC" = rgb(252, 97, 213, maxColorValue = 255),
  "Prog_Mk" = rgb(231, 107, 243, maxColorValue = 255),
  "Prog_RBC" = rgb(220, 113, 250, maxColorValue = 255),
  "Treg" = rgb(243, 123, 89, maxColorValue = 255)
)

plot(1:length(col_palette), pch = 16, col = col_palette, cex = 5)

