round4_config <- list(
 seed=20260905L,harmony_variable="orig.ident",resolutions=c(0.4,1.2),npcs=50L,nfeatures=3000L,
 input_patterns=c(Myeloid="Myeloid_round3_reclustered_preannotation.qs2",Fibroblast="Fibroblast_round3_reclustered_preannotation.qs2",Endothelial="Endothelial_round3_reclustered_preannotation.qs2",TNK="TNK_round3_reclustered_preannotation.qs2",Mural="Mural_round3_reclustered_preannotation.qs2",Schwann="Schwann_round3_extracted.qs2"),
 point_sizes=c(Fibroblast=.55,Mural=.70,TNK=.25)
)
round4_res_field<-function(x)paste0("round4_res_",format(x,nsmall=1,trim=TRUE))
round4_palette<-c("#E41A1C","#377EB8","#008837","#984EA3","#FF7F00","#A65628","#C51B7D","#1B9E77","#D95F02","#3F007D","#0072B2","#D73027","#238B45","#54278F","#8C510A","#01665E","#B2182B","#2166AC","#762A83","#A50F15","#00441B","#67000D","#08306B","#4D004B","#7F2704","#006D2C","#6A51A3","#CB181D","#08519C","#993404","#004C6D","#7A5195","#EF5675","#FFA600","#003F5C")
