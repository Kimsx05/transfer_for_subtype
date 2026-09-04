#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(SeuratObject)})
stopifnot(requireNamespace("qs2",quietly=TRUE))
root<-normalizePath(file.path(getwd(),"subtype"),mustWork=TRUE);r3<-file.path(root,"round3")
source(file.path(r3,"code","00_round3_config.R"))
dirs<-c("00_input_audit","01_extraction/Schwann","01_extraction/Mural","02_reclustering/Myeloid",
 "02_reclustering/Fibroblast","02_reclustering/Endothelial","02_reclustering/TNK","02_reclustering/Mural",
 "03_transition_audit","04_summary","logs")
invisible(lapply(file.path(r3,dirs),dir.create,recursive=TRUE,showWarnings=FALSE))

audit<-list();provenance<-list();objects<-list()
for(lineage in names(round3_config$inputs)){
 path<-file.path(root,round3_config$inputs[[lineage]]);obj<-qs2::qs_read(path);meta<-obj[[]]
 fld<-round3_config$working_fields[[lineage]]
 if(!fld%in%names(meta))stop("Missing frozen working field: ",lineage," / ",fld)
 if(!"orig.ident"%in%names(meta)||anyNA(meta$orig.ident))stop("orig.ident absent/missing: ",lineage)
 assays<-Assays(obj);layers<-unlist(lapply(assays,function(a)paste0(a,":",paste(Layers(obj[[a]],search=NA),collapse=";"))))
 audit[[lineage]]<-data.frame(lineage=lineage,object_path=path,n_cells=ncol(obj),assays=paste(assays,collapse=";"),
   default_assay=DefaultAssay(obj),layers=paste(layers,collapse=" | "),reductions=paste(Reductions(obj),collapse=";"),
   round2_resolution_field=fld,cluster_levels=paste(sort(unique(as.character(meta[[fld]]))),collapse=";"),
   metadata_fields=paste(names(meta),collapse=";"))
 provenance[[lineage]]<-data.frame(cell_barcode=rownames(meta),round2_compartment=lineage,
   round2_resolution=fld,round2_cluster=as.character(meta[[fld]]),round3_action="KEEP_MAIN_LINEAGE",
   round3_destination=lineage,stringsAsFactors=FALSE)
 objects[[lineage]]<-obj
}
write.csv(do.call(rbind,audit),file.path(r3,"00_input_audit","round2_input_object_audit.csv"),row.names=FALSE)

fib<-objects$Fibroblast;fm<-fib[[]];ff<-round3_config$working_fields[["Fibroblast"]]
schwann_cells<-rownames(fm)[as.character(fm[[ff]])=="15"]
mural_cells<-rownames(fm)[as.character(fm[[ff]])%in%c("3","4","13")]
if(!length(schwann_cells)||!length(mural_cells)||length(intersect(schwann_cells,mural_cells)))stop("Invalid Schwann/Mural partition")
provenance$Fibroblast$round3_action[provenance$Fibroblast$cell_barcode%in%schwann_cells]<-"EXTRACT_SCHWANN"
provenance$Fibroblast$round3_destination[provenance$Fibroblast$cell_barcode%in%schwann_cells]<-"Schwann"
provenance$Fibroblast$round3_action[provenance$Fibroblast$cell_barcode%in%mural_cells]<-"EXTRACT_MURAL"
provenance$Fibroblast$round3_destination[provenance$Fibroblast$cell_barcode%in%mural_cells]<-"Mural"

schwann<-subset(fib,cells=schwann_cells);mural<-subset(fib,cells=mural_cells)
schwann$round3_action<-"EXTRACT_SCHWANN";schwann$round3_destination<-"Schwann";schwann$round3_identity_confirmation<-"Schwann Cell"
mural$round3_action<-"EXTRACT_MURAL";mural$round3_destination<-"Mural"
mural$round2_preliminary_identity<-ifelse(as.character(mural[[ff,drop=TRUE]])=="4","Pericyte-like","Smooth Muscle / contractile mural-like")
qs2::qs_save(schwann,file.path(r3,"01_extraction","Schwann","Schwann_round3_extracted.qs2"),nthreads=16)
qs2::qs_save(mural,file.path(r3,"01_extraction","Mural","Mural_round3_extracted_raw.qs2"),nthreads=16)

manifest_fields<-function(obj,extra=list()){
 m<-obj[[]];wanted<-intersect(c("orig.ident","dataset","sample_id","patient_uid","stage","Type","T_stage","N_stage","M_stage","major_celltype"),names(m))
 out<-data.frame(cell_barcode=rownames(m),m[,wanted,drop=FALSE],check.names=FALSE)
 for(n in names(extra))out[[n]]<-extra[[n]]
 out
}
sm<-manifest_fields(schwann,list(round2_cluster=as.character(schwann[[ff,drop=TRUE]]),round2_resolution=ff))
write.csv(sm,file.path(r3,"01_extraction","Schwann","Schwann_cell_manifest.csv"),row.names=FALSE)
schwann_markers<-c("SOX10","PLP1","MPZ","PMP2","CDH19","S100B","ERBB3","GFRA3")
DefaultAssay(schwann)<-"RNA";present<-schwann_markers%in%rownames(schwann)
ma<-data.frame(gene=schwann_markers,present=present,stringsAsFactors=FALSE)
if(any(present)){
 dat<-FetchData(schwann,vars=schwann_markers[present],layer="data")
 ma$average_normalized_expression[match(colnames(dat),ma$gene)]<-colMeans(dat)
 ma$percent_expressed[match(colnames(dat),ma$gene)]<-100*colMeans(dat>0)
}
write.csv(ma,file.path(r3,"01_extraction","Schwann","Schwann_marker_audit.csv"),row.names=FALSE)
sms<-do.call(rbind,lapply(intersect(c("orig.ident","dataset","sample_id","patient_uid","stage","Type","major_celltype"),names(schwann[[]])),function(f){
 z<-as.data.frame(table(value=schwann[[f,drop=TRUE]],useNA="ifany"),stringsAsFactors=FALSE);names(z)[2]<-"n_cells";z$field<-f;z[,c("field","value","n_cells")]
}))
write.csv(sms,file.path(r3,"01_extraction","Schwann","Schwann_metadata_summary.csv"),row.names=FALSE)
mm<-manifest_fields(mural,list(round2_cluster=as.character(mural[[ff,drop=TRUE]]),round2_resolution=ff,
 round2_preliminary_identity=mural$round2_preliminary_identity))
write.csv(mm,file.path(r3,"01_extraction","Mural","Mural_cell_manifest.csv"),row.names=FALSE)

fib_retained<-setdiff(rownames(fm),c(schwann_cells,mural_cells))
part<-data.frame(round2_fibro_total=nrow(fm),schwann_extracted=length(schwann_cells),mural_extracted=length(mural_cells),
 round3_fibro_retained=length(fib_retained),accounting_ok=nrow(fm)==length(schwann_cells)+length(mural_cells)+length(fib_retained))
write.csv(part,file.path(r3,"01_extraction","Fibro_partition_audit.csv"),row.names=FALSE)
if(!part$accounting_ok)stop("Fibroblast partition accounting failed")
prov<-do.call(rbind,provenance)
if(nrow(prov)!=sum(vapply(objects,ncol,numeric(1)))||anyDuplicated(prov$cell_barcode))stop("Provenance incomplete/duplicated")
if(!all(prov$round3_action%in%c("KEEP_MAIN_LINEAGE","EXTRACT_SCHWANN","EXTRACT_MURAL")))stop("Invalid action")
write.csv(prov,file.path(r3,"04_summary","round3_cell_provenance_manifest.csv"),row.names=FALSE)
exsum<-data.frame(destination=c("Schwann","Mural","Fibroblast retained"),n_cells=c(length(schwann_cells),length(mural_cells),length(fib_retained)))
write.csv(exsum,file.path(r3,"04_summary","round3_extraction_summary.csv"),row.names=FALSE)
message("INPUT AUDIT AND EXTRACTION COMPLETE")
