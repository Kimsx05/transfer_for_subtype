round3_config <- list(
  seed = 20260904L,
  harmony_variable = "orig.ident",
  resolutions = c(0.4, 0.8, 1.2),
  npcs = 50L,
  nfeatures = 3000L,
  inputs = c(
    Myeloid = "results/Myeloid/Myeloid_round2_reclustered_preannotation.qs2",
    Fibroblast = "results/Fibroblast/Fibroblast_round2_reclustered_preannotation.qs2",
    Endothelial = "results/Endothelial/Endothelial_round2_reclustered_preannotation.qs2",
    TNK = "results/TNK/TNK_round2_reclustered_preannotation.qs2"
  ),
  working_fields = c(
    Myeloid = "round2_res_0.4",
    Fibroblast = "round2_res_0.8",
    Endothelial = "round2_res_0.6",
    TNK = "round2_res_0.8"
  ),
  point_sizes = c(Myeloid=.50, Fibroblast=.55, Endothelial=.60, TNK=.25, Mural=.70)
)

round3_res_field <- function(x) paste0("round3_res_", format(x, nsmall=1, trim=TRUE))
round3_palette <- c(
  "#E41A1C","#377EB8","#008837","#984EA3","#FF7F00","#A65628","#C51B7D","#1B9E77",
  "#D95F02","#3F007D","#0072B2","#D73027","#238B45","#54278F","#8C510A","#01665E",
  "#B2182B","#2166AC","#762A83","#A50F15","#00441B","#67000D","#08306B","#4D004B",
  "#7F2704","#006D2C","#6A51A3","#CB181D","#08519C","#993404","#004C6D","#7A5195",
  "#EF5675","#FFA600","#003F5C","#2F4B7C","#665191","#A05195","#D45087","#F95D6A"
)
