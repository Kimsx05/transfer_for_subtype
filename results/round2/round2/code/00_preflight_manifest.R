#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(SeuratObject)})
stopifnot(requireNamespace("qs",quietly=TRUE),requireNamespace("digest",quietly=TRUE))
wd<-normalizePath(getwd(),mustWork=TRUE); r2<-wd; r1<-normalizePath(file.path(wd,"..","round1"),mustWork=TRUE)
dir.create(file.path(r2,"audit"),recursive=TRUE,showWarnings=FALSE);dir.create(file.path(r2,"logs"),recursive=TRUE,showWarnings=FALSE)
spec<-list(
 B=list(md="B_round2_reclustering_plan_res0.6.md",exclude=c(12,13,14,16,17),resolution=.6,harmony="orig.ident; theta/default; lambda/default"),
 Endothelial=list(md="Endothelial_round2_reclustering_plan_res0.6.md",exclude=c(7,11,18,20,22),resolution=.6,harmony="orig.ident; theta/default; lambda/default"),
 Fibroblast=list(md="Fibroblast_round2_reclustering_plan_res1.0.md",exclude=c(10,11,12,19,21,23),resolution=1,harmony="orig.ident; theta/default; lambda/default"),
 Myeloid=list(md="Myeloid_round2_reclustering_plan_res0.8.md",exclude=c(5,11,19,23),resolution=.8,harmony="orig.ident; theta/default; lambda/default"),
 TNK=list(md="TNK_round2_harmony_sensitivity_res0.8.md",exclude=integer(),resolution=.8,harmony="A orig.ident theta2 lambda1; B orig.ident theta4 lambda1; C orig.ident theta6 lambda0.5; D dataset theta4 lambda1; max_iter_harmony=20"))
rows<-list(); issues<-character()
for(nm in names(spec)){
 s<-spec[[nm]]; md<-file.path(r2,s$md); q<-file.path(r1,nm,paste0(nm,"_round1_clustered.qs"))
 if(!file.exists(md))issues<-c(issues,paste("missing MD",md));if(!file.exists(q))issues<-c(issues,paste("missing QS",q))
 if(!file.exists(md)||!file.exists(q))next
 o<-qs::qread(q); req<-c("cluster_v2_round1","orig.ident","dataset"); miss<-setdiff(req,colnames(o[[]]));if(length(miss))issues<-c(issues,paste(nm,"missing metadata",paste(miss,collapse=",")))
 ly<-Layers(o[["RNA"]],search=NA);if(!all(c("counts","data")%in%ly))issues<-c(issues,paste(nm,"missing RNA counts/data"))
 cl<-sort(unique(as.integer(as.character(o$cluster_v2_round1)))); absent<-setdiff(s$exclude,cl);if(length(absent))issues<-c(issues,paste(nm,"exclusion clusters absent",paste(absent,collapse=",")))
 keep<-!as.character(o$cluster_v2_round1)%in%as.character(s$exclude)
 rows[[nm]]<-data.frame(lineage=nm,round1_qs=normalizePath(q),round1_cells=ncol(o),excluded_clusters=paste(s$exclude,collapse=","),excluded_cells=sum(!keep),round2_input_cells=sum(keep),resolution=s$resolution,harmony=s$harmony,marker_change=if(nm=="Fibroblast")"append only MyoFibro / SMC: ACTA2,TAGLN,MYH11" else "none",md=s$md,md5=digest::digest(file=md,algo="md5"),output_dir=file.path(r2,nm),stringsAsFactors=FALSE)
 rm(o);invisible(gc())
}
tab<-do.call(rbind,rows);write.csv(tab,file.path(r2,"audit","preflight_inputs.csv"),row.names=FALSE)
lines<-c("# Round2 execution manifest","",paste0("Generated UTC: ",format(Sys.time(),tz="UTC",usetz=TRUE)),"","## Frozen MD files","",paste0("- `",tab$md,"` (md5: `",tab$md5,"` )"),"","## Execution matrix","",paste0("- ",tab$lineage,": Round1 input `",tab$round1_qs,"`; input=",tab$round1_cells,"; exclusion=[",tab$excluded_clusters,"]; excluded cells=",tab$excluded_cells,"; Round2 input=",tab$round2_input_cells,"; resolution=",tab$resolution,"; Harmony=",tab$harmony,"; marker changes=",tab$marker_change,"; output=`",tab$output_dir,"`"),"","## Conflict and completeness audit","",if(length(issues))paste0("- CONFLICT/MISSING: ",issues) else "- No substantive conflicts or missing required inputs detected.","- Round1 objects are read-only inputs; Round1 cluster metadata will be retained.","- Round2 exclusions apply only to lineage working copies.","- Stopping point: after reclustering, visualization, markers, composition, mapping, and TNK sensitivity comparison; before annotation/cell2location.")
writeLines(lines,file.path(r2,"round2_execution_manifest.md"));writeLines(capture.output(sessionInfo()),file.path(r2,"audit","preflight_sessionInfo.txt"))
if(length(issues))stop(paste(issues,collapse=" | "))
cat("PREFLIGHT_OK\n")
