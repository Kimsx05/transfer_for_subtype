#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(SeuratObject)})
root<-normalizePath(file.path(getwd(),"subtype"),mustWork=TRUE);r3<-file.path(root,"round3")
source(file.path(r3,"code","00_round3_config.R"))
tracked<-list(Myeloid=c("8"),Endothelial=c("3","7"),TNK=c("1","15"),Fibroblast=c("8","11","12"))
short<-c(Myeloid="Myeloid",Endothelial="Endo",TNK="TNK",Fibroblast="Fibro")
split_rows<-list();transition_summary<-list();k<-0L
entropy<-function(p){p<-p[p>0];-sum(p*log2(p))}
for(lineage in names(tracked)){
 old<-qs2::qs_read(file.path(root,round3_config$inputs[[lineage]]));oldmeta<-old[[]];oldfield<-round3_config$working_fields[[lineage]]
 new<-qs2::qs_read(file.path(r3,"02_reclustering",lineage,paste0(lineage,"_round3_reclustered_preannotation.qs2")));newmeta<-new[[]]
 for(cl in tracked[[lineage]]){
   oldcells<-rownames(oldmeta)[as.character(oldmeta[[oldfield]])==cl];cells<-intersect(oldcells,rownames(newmeta))
   if(length(cells)!=length(oldcells))stop("Tracked cells missing: ",lineage," C",cl)
   rows<-do.call(rbind,lapply(round3_config$resolutions,function(rr){z<-as.data.frame(table(round3_cluster=newmeta[cells,round3_res_field(rr)]),stringsAsFactors=FALSE);names(z)[2]<-"n_cells";z<-z[z$n_cells>0,,drop=FALSE];z$round2_cluster=cl;z$round3_resolution=rr;z$fraction_of_round2_cluster=z$n_cells/length(cells);z[,c("round2_cluster","round3_resolution","round3_cluster","n_cells","fraction_of_round2_cluster")]}))
   write.csv(rows,file.path(r3,"03_transition_audit",paste0(short[[lineage]],"_R2C",cl,"_to_R3_transition.csv")),row.names=FALSE)
   counts<-sapply(round3_config$resolutions,function(rr)sum(rows$round3_resolution==rr))
   p<-sort(rows$fraction_of_round2_cluster[rows$round3_resolution==1.2],decreasing=TRUE)
   substantial<-sum(p>=.10);status<-if(length(p)&&p[1]>=.80)"STABLE_CLUSTER" else if(substantial>=2&&sum(head(p,3))>=.50)"NATURAL_SPLIT_SUPPORTED" else "DISTRIBUTED_REVIEW"
   k<-k+1L;split_rows[[k]]<-data.frame(lineage=lineage,round2_cluster=cl,round2_n_cells=length(cells),n_round3_clusters_at_0.4=counts[1],n_round3_clusters_at_0.8=counts[2],n_round3_clusters_at_1.2=counts[3],largest_fraction_at_1.2=if(length(p))p[1]else NA,second_largest_fraction_at_1.2=if(length(p)>1)p[2]else 0,split_entropy_at_1.2=entropy(p),split_status=status)
 }
 rm(old,new);gc()
}
splits<-do.call(rbind,split_rows);write.csv(splits,file.path(r3,"03_transition_audit","round3_tracked_cluster_split_summary.csv"),row.names=FALSE);write.csv(splits,file.path(r3,"04_summary","round3_transition_summary.csv"),row.names=FALSE)

objs<-c("Myeloid","Fibroblast","Endothelial","TNK","Mural");summaries<-list();valid<-list()
for(lineage in objs){base<-file.path(r3,"02_reclustering",lineage);s<-read.csv(file.path(base,"tables","round3_object_summary.csv"),check.names=FALSE);summaries[[lineage]]<-s;objfile<-file.path(base,paste0(lineage,"_round3_reclustered_preannotation.qs2"));obj<-qs2::qs_read(objfile);meta<-obj[[]];marker_ok<-sapply(round3_config$resolutions,function(rr)all(file.exists(file.path(base,"tables",paste0("Round3_markers_res",rr,c("_full.csv","_top100_raw.csv","_top100_annotation_friendly.csv","_top100_non_cellcycle.csv"))))));plot_ok<-sapply(round3_config$resolutions,function(rr)all(file.exists(file.path(base,"figures",paste0("Round3_",lineage,"_UMAP_res",rr,c(".pdf",".png"))),file.path(base,"figures",paste0("Round3_",lineage,"_canonical_DotPlot_res",rr,c(".pdf",".png"))))));valid[[lineage]]<-data.frame(object=lineage,n_cells=ncol(obj),object_file=objfile,reductions_ok=all(c("round3_pca","round3_harmony","round3_umap")%in%Reductions(obj)),resolution_fields_ok=all(round3_res_field(round3_config$resolutions)%in%names(meta)),harmony_orig_ident=all(obj$round3_harmony_variable=="orig.ident"),no_final_annotation_fields=!any(c("fine_celltype","final_celltype")%in%names(meta)),all_marker_files=all(marker_ok),all_plot_files=all(plot_ok));rm(obj);gc()}
write.csv(do.call(rbind,summaries),file.path(r3,"04_summary","round3_all_objects_summary.csv"),row.names=FALSE)
validation<-do.call(rbind,valid);write.csv(validation,file.path(r3,"04_summary","round3_marker_output_validation.csv"),row.names=FALSE)
if(!all(unlist(validation[,c("reductions_ok","resolution_fields_ok","harmony_orig_ident","no_final_annotation_fields","all_marker_files","all_plot_files")])))stop("Round3 validation failed")
prov<-read.csv(file.path(r3,"04_summary","round3_cell_provenance_manifest.csv"),stringsAsFactors=FALSE)
if(anyDuplicated(prov$cell_barcode)||!all(prov$round3_action%in%c("KEEP_MAIN_LINEAGE","EXTRACT_SCHWANN","EXTRACT_MURAL")))stop("Provenance validation failed")
cat("ROUND3 COMPLETE\n\nSCHWANN EXTRACTED\nMURAL EXTRACTED AND RECLUSTERED\nNO LOCAL RESCUE PERFORMED\nNO AMBIGUOUS MIXED CLUSTERS DELETED\nMAIN LINEAGES RECLUSTERED AT 0.4 / 0.8 / 1.2\nROUND2 SUSPICIOUS CLUSTERS TRACKED INTO ROUND3\n\nSTOPPED BEFORE FINAL RESOLUTION SELECTION AND ANNOTATION\n")
