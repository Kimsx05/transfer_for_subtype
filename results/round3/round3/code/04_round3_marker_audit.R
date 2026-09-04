#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(ggplot2)})
if(requireNamespace("future",quietly=TRUE))future::plan("sequential")
options(future.globals.maxSize=200*1024^3)
root<-normalizePath(file.path(getwd(),"subtype"),mustWork=TRUE);r3<-file.path(root,"round3")
source(file.path(r3,"code","00_round3_config.R"));source(file.path(r3,"code","01_round3_marker_panels.R"))
lineage<-Sys.getenv("ROUND3_OBJECT",unset="");if(!lineage%in%names(round3_marker_panels))stop("Set ROUND3_OBJECT")
out<-file.path(r3,"02_reclustering",lineage);tab<-file.path(out,"tables");fig<-file.path(out,"figures")
con<-file(file.path(r3,"logs",paste0("Round3_",lineage,"_marker_audit.log")),"wt");sink(con,type="output",split=TRUE);sink(con,type="message")
message("START Round3 marker audit ",lineage,"; RNA normalized expression only")
object<-qs2::qs_read(file.path(out,paste0(lineage,"_round3_reclustered_preannotation.qs2")))
DefaultAssay(object)<-"RNA";layers<-Layers(object[["RNA"]],search=NA);data_layers<-grep("^data(\\.|$)",layers,value=TRUE)
if(length(data_layers)>1L){message("JoinLayers in marker-only memory copy");object<-JoinLayers(object,assay="RNA")}else message("JoinLayers not required; layers = ",paste(layers,collapse=","))
cycle<-c("MKI67","TOP2A","UBE2C","BIRC5","CENPF","RRM2","PCLAF","CDC20","TYMS","CDK1","CCNB1","CCNB2","CCNA2","CDCA3","CDCA8","CDCA5","CDCA2","TK1","NUSAP1","ASPM","TPX2","AURKA","AURKB","KIF11","KIF20A","KIF2C","KIF4A","KIFC1","CENPA","CENPE","CENPU","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","PCNA","HMGB2","STMN1","TUBA1B")
top_n<-function(x,n=100,mode=c("raw","friendly","noncycle")){mode<-match.arg(mode);z<-x;if(mode=="friendly")z<-z[!grepl("^(MT-|RPL|RPS|ENSG)",z$gene)&z$gene!="MALAT1",,drop=FALSE];if(mode=="noncycle")z<-z[!z$gene%in%cycle&!grepl("^HIST[1-4]",z$gene),,drop=FALSE];fc<-intersect(c("avg_log2FC","avg_logFC"),names(z));if(!length(fc))stop("FC column absent");z<-z[order(as.character(z$cluster),-z[[fc[1]]],z$p_val_adj),,drop=FALSE];do.call(rbind,lapply(split(z,z$cluster),head,n=n))}
summary<-list()
for(rr in round3_config$resolutions){f<-round3_res_field(rr);if(!f%in%names(object[[]]))stop("Missing ",f);Idents(object)<-object[[f,drop=TRUE]];message("FindAllMarkers ",lineage," ",f);m<-FindAllMarkers(object,assay="RNA",slot="data",only.pos=TRUE,min.pct=.10,logfc.threshold=.25,test.use="wilcox",densify=FALSE,verbose=TRUE);raw<-top_n(m,100,"raw");friendly<-top_n(m,100,"friendly");noncycle<-top_n(m,100,"noncycle");stem<-file.path(tab,paste0("Round3_markers_res",rr));write.csv(m,paste0(stem,"_full.csv"),row.names=FALSE);write.csv(raw,paste0(stem,"_top100_raw.csv"),row.names=FALSE);write.csv(friendly,paste0(stem,"_top100_annotation_friendly.csv"),row.names=FALSE);write.csv(noncycle,paste0(stem,"_top100_non_cellcycle.csv"),row.names=FALSE);summary[[as.character(rr)]]<-data.frame(object=lineage,resolution=rr,object_clusters=length(unique(object[[f,drop=TRUE]])),marker_clusters=length(unique(m$cluster)),full_rows=nrow(m),top100_raw_rows=nrow(raw),top100_friendly_rows=nrow(friendly),top100_noncycle_rows=nrow(noncycle));rm(m,raw,friendly,noncycle);invisible(gc())}
write.csv(do.call(rbind,summary),file.path(tab,"Round3_marker_output_summary.csv"),row.names=FALSE)
panels<-round3_marker_panels[[lineage]];panel<-do.call(rbind,lapply(names(panels),function(p)data.frame(program=p,gene=panels[[p]],stringsAsFactors=FALSE)));panel$present<-panel$gene%in%rownames(object);write.csv(panel,file.path(tab,"Round3_canonical_marker_panel_audit.csv"),row.names=FALSE);missing<-panel[!panel$present,,drop=FALSE];write.csv(missing,file.path(tab,"Round3_missing_canonical_markers.csv"),row.names=FALSE)
features<-unique(panel$gene[panel$present])
for(rr in round3_config$resolutions){f<-round3_res_field(rr);Idents(object)<-object[[f,drop=TRUE]];p<-DotPlot(object,features=features,assay="RNA",group.by=f,scale=TRUE,dot.scale=6)+RotatedAxis()+labs(x="marker",y="cluster",title=paste0("res = ",rr))+theme_classic(base_size=10)+theme(plot.title=element_text(hjust=.5),axis.text.x=element_text(size=7));stem<-file.path(fig,paste0("Round3_",lineage,"_canonical_DotPlot_res",rr));w<-max(14,.25*length(features));ggsave(paste0(stem,".pdf"),p,width=w,height=8,limitsize=FALSE);ggsave(paste0(stem,".png"),p,width=w,height=8,dpi=300,bg="white",limitsize=FALSE)}
message("END Round3 marker audit ",lineage,"; NO ANNOTATION")
sink(type="message");sink(type="output");close(con)
