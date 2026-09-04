#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(ggplot2);library(patchwork)})
stopifnot(requireNamespace("qs2",quietly=TRUE),requireNamespace("harmony",quietly=TRUE),requireNamespace("ggrepel",quietly=TRUE))
if(requireNamespace("future",quietly=TRUE))future::plan("sequential")
options(future.globals.maxSize=200*1024^3)
root<-normalizePath(file.path(getwd(),"subtype"),mustWork=TRUE);r3<-file.path(root,"round3")
source(file.path(r3,"code","00_round3_config.R"))
lineage<-Sys.getenv("ROUND3_OBJECT",unset="")
if(!lineage%in%c("Myeloid","Fibroblast","Endothelial","TNK","Mural"))stop("Set ROUND3_OBJECT")
set.seed(round3_config$seed);resolutions<-round3_config$resolutions
out<-file.path(r3,"02_reclustering",lineage);tab<-file.path(out,"tables");fig<-file.path(out,"figures")
dir.create(tab,recursive=TRUE,showWarnings=FALSE);dir.create(fig,recursive=TRUE,showWarnings=FALSE)
logcon<-file(file.path(r3,"logs",paste0("Round3_",lineage,"_reclustering.log")),"wt")
sink(logcon,type="output",split=TRUE);sink(logcon,type="message")
message("START Round3 ",lineage,"; Harmony variable = orig.ident; resolutions = 0.4/0.8/1.2")

if(lineage=="Mural"){
 input<-file.path(r3,"01_extraction","Mural","Mural_round3_extracted_raw.qs2");object<-qs2::qs_read(input)
}else{
 input<-file.path(root,round3_config$inputs[[lineage]]);object<-qs2::qs_read(input)
 if(lineage=="Fibroblast"){
   prov<-read.csv(file.path(r3,"04_summary","round3_cell_provenance_manifest.csv"),stringsAsFactors=FALSE)
   keep<-prov$cell_barcode[prov$round2_compartment=="Fibroblast"&prov$round3_action=="KEEP_MAIN_LINEAGE"]
   object<-subset(object,cells=keep)
 }
 object$round3_action<-"KEEP_MAIN_LINEAGE";object$round3_destination<-lineage
}
n_input_cells<-if(lineage=="Fibroblast") ncol(object)+
  nrow(read.csv(file.path(r3,"01_extraction","Schwann","Schwann_cell_manifest.csv")))+
  nrow(read.csv(file.path(r3,"01_extraction","Mural","Mural_cell_manifest.csv"))) else ncol(object)
if(!"orig.ident"%in%names(object[[]])||anyNA(object$orig.ident))stop("orig.ident missing")
object$round3_source_object<-input;object$round3_harmony_variable<-"orig.ident"
object@reductions<-list();object@graphs<-list();object@neighbors<-list();object@commands<-list()
DefaultAssay(object)<-"RNA";layers<-Layers(object[["RNA"]],search=NA);counts<-grep("^counts(\\.|$)",layers,value=TRUE)
if(!length(counts))stop("RNA counts layer absent")
if(length(counts)>1L){message("JoinLayers on Round3 object copy only; counts layers = ",length(counts));object<-JoinLayers(object,assay="RNA")
}else message("JoinLayers not required; RNA layers before normalization = ",paste(layers,collapse=","))
object<-NormalizeData(object,assay="RNA",normalization.method="LogNormalize",scale.factor=10000,verbose=TRUE)
object<-FindVariableFeatures(object,assay="RNA",selection.method="vst",nfeatures=min(round3_config$nfeatures,nrow(object)),verbose=TRUE)
object<-ScaleData(object,assay="RNA",features=VariableFeatures(object),verbose=TRUE)
npcs<-min(round3_config$npcs,length(VariableFeatures(object))-1L,ncol(object)-1L)
if(npcs<20)stop("Insufficient dimensions/cells")
object<-RunPCA(object,assay="RNA",features=VariableFeatures(object),npcs=npcs,reduction.name="round3_pca",
 reduction.key="R3PC_",seed.use=round3_config$seed,verbose=TRUE)
stdev<-Stdev(object[["round3_pca"]]);pct<-100*stdev^2/sum(stdev^2);cum<-cumsum(pct)
candidates<-which(cum>=80&pct<5);raw<-if(length(candidates))candidates[[1]]else length(stdev)
lower<-if(ncol(object)<2000)15L else 20L;upper<-if(ncol(object)<5000)35L else 40L
ndims<-as.integer(max(lower,min(upper,raw,length(stdev))))
pca<-data.frame(PC=seq_along(stdev),stdev=stdev,variance_percent=pct,cumulative_variance_percent=cum,
 selected_for_round3=seq_along(stdev)<=ndims)
write.csv(pca,file.path(tab,paste0("Round3_",lineage,"_pca_variance_audit.csv")),row.names=FALSE)
pe<-ggplot(pca,aes(PC,variance_percent))+geom_point(size=1)+geom_line()+geom_vline(xintercept=ndims,linetype=2,colour="red")+
 theme_classic()+labs(title=paste0("selected dims = 1:",ndims),y="Variance (%)")
ggsave(file.path(fig,paste0("Round3_",lineage,"_PCA_elbow.pdf")),pe,width=8,height=6)
ggsave(file.path(fig,paste0("Round3_",lineage,"_PCA_elbow.png")),pe,width=8,height=6,dpi=300,bg="white")
message("PCA dimensions = ",length(stdev),"; selected dims = 1:",ndims,"; cells = ",ncol(object))
object<-harmony::RunHarmony(object,group.by.vars="orig.ident",reduction.use="round3_pca",dims.use=seq_len(ndims),
 reduction.save="round3_harmony",verbose=TRUE)
object<-FindNeighbors(object,reduction="round3_harmony",dims=seq_len(ndims),graph.name=c("round3_nn","round3_snn"),verbose=TRUE)
for(rr in resolutions)object<-FindClusters(object,graph.name="round3_snn",resolution=rr,cluster.name=round3_res_field(rr),
 algorithm=1,random.seed=round3_config$seed,group.singletons=FALSE,verbose=TRUE)
object<-RunUMAP(object,reduction="round3_harmony",dims=seq_len(ndims),reduction.name="round3_umap",reduction.key="R3UMAP_",
 seed.use=round3_config$seed,verbose=TRUE)
object$round3_clustering_dims<-ndims;object$round3_annotation_status<-"preannotation"
outfile<-file.path(out,paste0(lineage,"_round3_reclustered_preannotation.qs2"));qs2::qs_save(object,outfile,nthreads=16)

sizes<-do.call(rbind,lapply(resolutions,function(rr){z<-as.data.frame(table(cluster=object[[round3_res_field(rr),drop=TRUE]]),stringsAsFactors=FALSE);names(z)[2]<-"n_cells";z$resolution<-rr;z$fraction<-z$n_cells/sum(z$n_cells);z[,c("resolution","cluster","n_cells","fraction")]}))
write.csv(sizes,file.path(tab,"round3_resolution_cluster_sizes.csv"),row.names=FALSE)
sums<-do.call(rbind,lapply(split(sizes,sizes$resolution),function(z)data.frame(resolution=unique(z$resolution),n_clusters=nrow(z),min_cluster_n=min(z$n_cells),median_cluster_n=median(z$n_cells),max_cluster_n=max(z$n_cells))))
write.csv(sums,file.path(tab,"round3_resolution_summary.csv"),row.names=FALSE)
trans<-do.call(rbind,lapply(1:2,function(i){a<-resolutions[i];b<-resolutions[i+1];z<-as.data.frame(table(parent_cluster=object[[round3_res_field(a),drop=TRUE]],child_cluster=object[[round3_res_field(b),drop=TRUE]]),stringsAsFactors=FALSE);names(z)[3]<-"n_cells";z<-z[z$n_cells>0,,drop=FALSE];z$from_resolution<-a;z$to_resolution<-b;z$fraction_of_parent<-ave(z$n_cells,z$parent_cluster,FUN=function(x)x/sum(x));z[,c("from_resolution","to_resolution","parent_cluster","child_cluster","n_cells","fraction_of_parent")]}))
write.csv(trans,file.path(tab,"round3_resolution_transition_table.csv"),row.names=FALSE)
for(i in 1:2){a<-resolutions[i];b<-resolutions[i+1];z<-as.data.frame.matrix(table(object[[round3_res_field(a),drop=TRUE]],object[[round3_res_field(b),drop=TRUE]]));write.csv(cbind(parent_cluster=rownames(z),z),file.path(tab,paste0("round3_transition_res",a,"_to_res",b,".csv")),row.names=FALSE)}

make_plot<-function(field,title,label=FALSE){emb<-Embeddings(object[["round3_umap"]]);d<-data.frame(x=emb[,1],y=emb[,2],group=as.character(object[[]][rownames(emb),field]));d$group[is.na(d$group)]<-"NA";d$group<-factor(d$group,levels=sort(unique(d$group)));cols<-setNames(rep(round3_palette,length.out=nlevels(d$group)),levels(d$group));p<-ggplot(d,aes(x,y,colour=group))+geom_point(size=round3_config$point_sizes[[lineage]],alpha=1,stroke=0)+scale_colour_manual(values=cols,drop=FALSE)+labs(x="UMAP_1",y="UMAP_2",colour=field,title=title)+theme_classic(base_size=11)+theme(plot.title=element_text(hjust=.5),legend.key.height=grid::unit(.35,"cm"));if(label){centers<-aggregate(cbind(x,y)~group,d,median);p<-p+ggrepel::geom_label_repel(data=centers,aes(label=group),colour="black",fill="white",size=3.3,fontface="bold",linewidth=.2,seed=round3_config$seed,max.overlaps=Inf,show.legend=FALSE)};p}
save_pair<-function(p,stem,w=10,h=8){ggsave(paste0(stem,".pdf"),p,width=w,height=h,limitsize=FALSE);ggsave(paste0(stem,".png"),p,width=w,height=h,dpi=300,bg="white",limitsize=FALSE)}
plots<-lapply(resolutions,function(rr)make_plot(round3_res_field(rr),paste0("res = ",rr),TRUE))
for(i in seq_along(resolutions))save_pair(plots[[i]],file.path(fig,paste0("Round3_",lineage,"_UMAP_res",resolutions[i])))
comparison<-plots[[1]]+plots[[2]]+plots[[3]]+plot_layout(ncol=3,guides="collect")
save_pair(comparison,file.path(fig,paste0("Round3_",lineage,"_UMAP_resolution_comparison")),21,7)
audit_fields<-intersect(c("dataset","orig.ident","stage","Type","T_stage","N_stage","M_stage"),names(object[[]]))
for(f in audit_fields)save_pair(make_plot(f,f,FALSE),file.path(fig,paste0("Round3_",lineage,"_UMAP_by_",f)))
composition<-list();enrich<-list();k<-0L
for(rr in resolutions)for(f in audit_fields){z<-as.data.frame(table(cluster=object[[round3_res_field(rr),drop=TRUE]],group=object[[f,drop=TRUE]],useNA="ifany"),stringsAsFactors=FALSE);names(z)[3]<-"n_cells";z$resolution<-rr;z$field<-f;z$within_cluster_fraction<-ave(z$n_cells,z$cluster,FUN=function(x)x/sum(x));k<-k+1L;composition[[k]]<-z;if(f=="dataset"){tot<-tapply(z$n_cells,z$group,sum);z$lineage_wide_dataset_fraction<-as.numeric(tot[as.character(z$group)])/sum(z$n_cells);z$enrichment<-z$within_cluster_fraction/z$lineage_wide_dataset_fraction;z$high_enrichment_low_count_flag<-z$enrichment>=3&z$n_cells<50;enrich[[as.character(rr)]]<-z}}
write.csv(do.call(rbind,composition),file.path(tab,"round3_cluster_metadata_composition.csv"),row.names=FALSE)
if(length(enrich))write.csv(do.call(rbind,enrich),file.path(tab,"round3_dataset_enrichment.csv"),row.names=FALSE)
summary<-data.frame(object=lineage,round2_n_cells=n_input_cells,n_extracted=if(lineage=="Fibroblast")n_input_cells-ncol(object)else 0,round3_n_cells=ncol(object),Harmony_variable="orig.ident",PCA_dims=length(stdev),neighbor_dims=ndims,n_clusters_res0.4=sums$n_clusters[sums$resolution==.4],n_clusters_res0.8=sums$n_clusters[sums$resolution==.8],n_clusters_res1.2=sums$n_clusters[sums$resolution==1.2],smallest_cluster_res0.4=sums$min_cluster_n[sums$resolution==.4],smallest_cluster_res0.8=sums$min_cluster_n[sums$resolution==.8],smallest_cluster_res1.2=sums$min_cluster_n[sums$resolution==1.2])
write.csv(summary,file.path(tab,"round3_object_summary.csv"),row.names=FALSE)
message("END Round3 ",lineage,"; NO LOCAL RESCUE; NO ANNOTATION")
sink(type="message");sink(type="output");close(logcon)
