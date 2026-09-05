#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(SeuratObject)})
stopifnot(requireNamespace("qs2",quietly=TRUE))
root<-normalizePath(file.path(getwd(),"subtype"),mustWork=TRUE);r4<-file.path(root,"round4_final");source(file.path(r4,"code","00_round4_config.R"))
dirs<-c("00_input_audit","01_lineage_correction/Fibro_to_Mural","01_lineage_correction/Mural_to_Fibro","01_lineage_correction/TNK_excluded","02_reclustering/Fibroblast","02_reclustering/Mural","02_reclustering/TNK","03_transition_audit","04_final_annotation_support","05_summary","logs")
invisible(lapply(file.path(r4,dirs),dir.create,recursive=TRUE,showWarnings=FALSE))
allqs<-list.files(file.path(root,"round3"),pattern="[.]qs2$",recursive=TRUE,full.names=TRUE)
paths<-sapply(round4_config$input_patterns,function(p){x<-allqs[basename(allqs)==p];if(length(x)!=1)stop("Expected exactly one input for ",p,"; found ",length(x));normalizePath(x)})
objects<-lapply(paths,qs2::qs_read);audit<-list()
for(n in names(objects)){o<-objects[[n]];m<-o[[]];r04<-grep("^round3_res_0[.]4$",names(m),value=TRUE);r12<-grep("^round3_res_1[.]2$",names(m),value=TRUE);if(n!="Schwann"&&(length(r04)!=1||length(r12)!=1))stop("Round3 fields not unique for ",n);ass<-Assays(o);lay<-unlist(lapply(ass,function(a)paste0(a,":",paste(Layers(o[[a]],search=NA),collapse=";"))));audit[[n]]<-data.frame(object=n,object_path=paths[[n]],n_cells=ncol(o),default_assay=DefaultAssay(o),assays=paste(ass,collapse=";"),layers=paste(lay,collapse=" | "),reductions=paste(Reductions(o),collapse=";"),res0.4_field=if(length(r04))r04 else NA,res1.2_field=if(length(r12))r12 else NA,cluster_levels=if(length(r12))paste(sort(unique(as.character(m[[r12]]))),collapse=";")else NA,metadata_fields=paste(names(m),collapse=";"))}
write.csv(do.call(rbind,audit),file.path(r4,"00_input_audit","round4_input_object_audit.csv"),row.names=FALSE)
f<-objects$Fibroblast;mu<-objects$Mural;t<-objects$TNK;ff<-"round3_res_1.2";mf<-"round3_res_1.2";tf<-"round3_res_1.2"
for(z in c("11","19","9"))if(!z%in%as.character(unique(f[[ff,drop=TRUE]])))stop("Missing Fibro C",z)
for(z in c("4","15"))if(!z%in%as.character(unique(mu[[mf,drop=TRUE]])))stop("Missing Mural C",z)
for(z in c("5","16"))if(!z%in%as.character(unique(t[[tf,drop=TRUE]])))stop("Missing TNK C",z)
cells_by<-function(o,field,cl)Cells(o)[as.character(o[[field,drop=TRUE]])==cl]
fc11<-cells_by(f,ff,"11");fc19<-cells_by(f,ff,"19");mc4<-cells_by(mu,mf,"4");mc15<-cells_by(mu,mf,"15");tc5<-cells_by(t,tf,"5");tc16<-cells_by(t,tf,"16");single<-cells_by(t,tf,"singleton")
subsave<-function(o,cells,path,action,dest){z<-subset(o,cells=cells);z$round4_action<-action;z$round4_destination<-dest;qs2::qs_save(z,path,nthreads=16);z}
f11<-subsave(f,fc11,file.path(r4,"01_lineage_correction","Fibro_to_Mural","Fibro_C11_to_Mural.qs2"),"MOVE_FIBRO_TO_MURAL","Mural")
f19<-subsave(f,fc19,file.path(r4,"01_lineage_correction","Fibro_to_Mural","Fibro_C19_to_Mural.qs2"),"MOVE_FIBRO_TO_MURAL","Mural")
m4<-subsave(mu,mc4,file.path(r4,"01_lineage_correction","Mural_to_Fibro","Mural_C4_to_Fibro.qs2"),"MOVE_MURAL_TO_FIBRO","Fibroblast")
m15<-subsave(mu,mc15,file.path(r4,"01_lineage_correction","Mural_to_Fibro","Mural_C15_to_Fibro.qs2"),"MOVE_MURAL_TO_FIBRO","Fibroblast")
e5<-subsave(t,tc5,file.path(r4,"01_lineage_correction","TNK_excluded","TNK_C5_excluded_from_reference.qs2"),"REMOVE_TNK_MIXED","EXCLUDED")
e16<-subsave(t,tc16,file.path(r4,"01_lineage_correction","TNK_excluded","TNK_C16_excluded_from_reference.qs2"),"REMOVE_TNK_MIXED","EXCLUDED")
if(length(single))es<-subsave(t,single,file.path(r4,"01_lineage_correction","TNK_excluded","TNK_singleton_excluded.qs2"),"REMOVE_SINGLETON","EXCLUDED")
manifest<-function(o,cells,source,field,action,dest,reason=NULL){m<-o[[]];ix<-match(cells,rownames(m));wanted<-intersect(c("orig.ident","dataset","sample_id","patient_uid","stage","Type","T_stage","N_stage","M_stage"),names(m));z<-data.frame(cell_barcode=cells,round3_source_object=source,round3_res1.2_cluster=as.character(m[ix,field]),round4_action=action,round4_destination=dest,m[ix,wanted,drop=FALSE],check.names=FALSE);if(!is.null(reason))z$exclusion_reason<-reason;z}
fm<-rbind(manifest(f,fc11,paths[["Fibroblast"]],ff,"MOVE_FIBRO_TO_MURAL","Mural"),manifest(f,fc19,paths[["Fibroblast"]],ff,"MOVE_FIBRO_TO_MURAL","Mural"));write.csv(fm,file.path(r4,"01_lineage_correction","Fibro_to_Mural","Fibro_to_Mural_manifest.csv"),row.names=FALSE)
mm<-rbind(manifest(mu,mc4,paths[["Mural"]],mf,"MOVE_MURAL_TO_FIBRO","Fibroblast"),manifest(mu,mc15,paths[["Mural"]],mf,"MOVE_MURAL_TO_FIBRO","Fibroblast"));write.csv(mm,file.path(r4,"01_lineage_correction","Mural_to_Fibro","Mural_to_Fibro_manifest.csv"),row.names=FALSE)
tm<-rbind(manifest(t,tc5,paths[["TNK"]],tf,"REMOVE_TNK_MIXED","EXCLUDED","EPITHELIAL_HEAVY_MIXED_TNK"),manifest(t,tc16,paths[["TNK"]],tf,"REMOVE_TNK_MIXED","EXCLUDED","STABLE_MULTI_LINEAGE_MIXED_POPULATION"),if(length(single))manifest(t,single,paths[["TNK"]],tf,"REMOVE_SINGLETON","EXCLUDED","GRAPH_SINGLETON"));write.csv(tm,file.path(r4,"01_lineage_correction","TNK_excluded","TNK_exclusion_manifest.csv"),row.names=FALSE)
fkeep<-setdiff(Cells(f),c(fc11,fc19));mkeep<-setdiff(Cells(mu),c(mc4,mc15));tkeep<-setdiff(Cells(t),c(tc5,tc16,single))
fbase<-subset(f,cells=fkeep);madd<-merge(m4,y=m15,merge.data=TRUE);finput<-merge(fbase,y=madd,merge.data=TRUE)
mubase<-subset(mu,cells=mkeep);fadd<-merge(f11,y=f19,merge.data=TRUE);muinput<-merge(mubase,y=fadd,merge.data=TRUE)
tinput<-subset(t,cells=tkeep)
for(pair in list(c("Fibroblast",ncol(finput)),c("Mural",ncol(muinput)),c("TNK",ncol(tinput))))message(pair[1]," input cells = ",pair[2])
if(anyDuplicated(Cells(finput))||anyDuplicated(Cells(muinput))||anyDuplicated(Cells(tinput)))stop("Duplicate barcode in Round4 input")
if(!"9"%in%as.character(f[[ff,drop=TRUE]][match(Cells(f)[as.character(f[[ff,drop=TRUE]])=="9"],Cells(f))]))stop("Fibro C9 check failed")
qs2::qs_save(finput,file.path(r4,"02_reclustering","Fibroblast","Fibroblast_round4_input.qs2"),nthreads=16);qs2::qs_save(muinput,file.path(r4,"02_reclustering","Mural","Mural_round4_input.qs2"),nthreads=16);qs2::qs_save(tinput,file.path(r4,"02_reclustering","TNK","TNK_round4_input.qs2"),nthreads=16)
prov<-list();makeprov<-function(name,o,source,field){m<-o[[]];data.frame(cell_barcode=Cells(o),round3_source=source,round3_res1.2_cluster=if(field%in%names(m))as.character(m[[field]])else NA,round4_action="KEEP",round4_destination=name,include_in_final_reference=TRUE,stringsAsFactors=FALSE)}
prov$Fibroblast<-makeprov("Fibroblast",f,paths[["Fibroblast"]],ff);prov$Mural<-makeprov("Mural",mu,paths[["Mural"]],mf);prov$TNK<-makeprov("TNK",t,paths[["TNK"]],tf)
setact<-function(z,cells,action,dest,include){ix<-z$cell_barcode%in%cells;z$round4_action[ix]<-action;z$round4_destination[ix]<-dest;z$include_in_final_reference[ix]<-include;z}
prov$Fibroblast<-setact(prov$Fibroblast,c(fc11,fc19),"MOVE_FIBRO_TO_MURAL","Mural",TRUE);prov$Mural<-setact(prov$Mural,c(mc4,mc15),"MOVE_MURAL_TO_FIBRO","Fibroblast",TRUE);prov$TNK<-setact(prov$TNK,c(tc5,tc16),"REMOVE_TNK_MIXED","EXCLUDED",FALSE);prov$TNK<-setact(prov$TNK,single,"REMOVE_SINGLETON","EXCLUDED",FALSE)
for(n in c("Myeloid","Endothelial","Schwann")){o<-objects[[n]];m<-o[[]];field<-if("round3_res_1.2"%in%names(m))"round3_res_1.2"else "";z<-makeprov(n,o,paths[[n]],field);z$round4_action<-"FROZEN_FROM_ROUND3";prov[[n]]<-z}
p<-do.call(rbind,prov);if(anyDuplicated(p$cell_barcode))stop("Duplicate provenance barcode");allowed<-c("KEEP","MOVE_FIBRO_TO_MURAL","MOVE_MURAL_TO_FIBRO","REMOVE_TNK_MIXED","REMOVE_SINGLETON","FROZEN_FROM_ROUND3");if(!all(p$round4_action%in%allowed))stop("Invalid provenance action");write.csv(p,file.path(r4,"05_summary","final_reference_cell_provenance.csv"),row.names=FALSE)
frozen<-data.frame(object=c("Myeloid","Endothelial","Schwann"),object_path=paths[c("Myeloid","Endothelial","Schwann")],status="FROZEN_FROM_ROUND3");write.csv(frozen,file.path(r4,"05_summary","frozen_object_manifest.csv"),row.names=FALSE)
acct<-data.frame(lineage=c("Fibroblast","Mural","TNK"),round3_cells=c(ncol(f),ncol(mu),ncol(t)),moved_in=c(length(mc4)+length(mc15),length(fc11)+length(fc19),0),moved_out=c(length(fc11)+length(fc19),length(mc4)+length(mc15),0),excluded=c(0,0,length(tc5)+length(tc16)+length(single)),round4_input_cells=c(ncol(finput),ncol(muinput),ncol(tinput)));acct$accounting_ok<-acct$round3_cells+acct$moved_in-acct$moved_out-acct$excluded==acct$round4_input_cells;write.csv(acct,file.path(r4,"05_summary","round4_lineage_correction_accounting.csv"),row.names=FALSE);if(!all(acct$accounting_ok))stop("Accounting failed")
message("ROUND4 INPUT AUDIT AND LINEAGE CORRECTION COMPLETE")
