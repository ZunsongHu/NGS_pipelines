options(future.globals.maxSize = 2*1024 * 1024^2) #10GB

# gc(rm(list=ls()))
options(warn=-1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(SummarizedExperiment)
  library(dplyr)
})
options(warn=0)


source("/net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/singleCell/scripts/R/singleR/singleR.R")

# functions   ------------------------------------------------------------------------------------- ---- 
get_barcode_df=function(object_in){
  data.frame(barcode=row.names(object_in@meta.data),stringsAsFactors = F)
}

dimPlot_fromDF=function(objIn = obj_QC,df = singleR_out_$df_singleR_out,
                        var = "labels",
                        label = paste0(sample,".",label,".Celllevel"),
                        cols_in = col_in,
                        dir_out = dir_out,w=12,h=9,
                        label_figure=T,repel_figure=T
){
  
  if(!all(colnames(obj_QC)==df$barcode)){stop("barcodes of obejct and df for var are not matched")}
  
  objIn[[label]]= unname(unlist(df[var]))
  
  DimPlot(objIn,group.by = label,label=label_figure,repel = repel_figure,
          cols = cols_in[names(cols_in) %in% unlist(objIn@meta.data[label])])
  ggsave(paste0(dir_out,label,".png"),width=w,height = h)
  ggsave(paste0(dir_out,label,".pdf"),width=w,height = h)
  objIn
}

# parameters temp   ------------------------------------------------------------------------------------- ---- 
# setwd("/net/nfs-irwrsrchnas01/labs/zgu_grp/Group/Grp_XiaojingYan/")
# 
# rds_in=     'analysis/obj/obj_seurat.SCT.TenSamples.rds'
# dir_out=    'analysis/singleR/sct/'
# sample=     "TenSamples"
# label=      "1M_20celltypes"
# ref_rds=    "/net/nfs-irwrsrchnas01/labs/zgu_grp/DryLab/refData/Celltype/human/boneMarrow/ref_1M_20Celltypes.rds"
# anno_label= "celltype"
# col_file=   "nan"
# transform=  "SCT"

# parameters  ------------------------------------------------------------------------------------- ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=args[1]
dir_out=args[2]
sample=args[3]
label=args[4]
ref_rds=args[5]
anno_label=args[6]
col_file=args[7]
transform=args[8]

# other parameters ------------------------------------------------------------------------------------- ----
do_subset=F;n_max=500
do_clustering=F

#create output folder  ------------------------------------------------------------------------------------- ----
if(!dir.exists(dir_out)){dir.create(dir_out,recursive = T)}

#read ref  ------------------------------------------------------------------------------------- ----
cat("Reading ref ...\n")
ref_=readRDS(ref_rds)
anno_levels=unique(unlist(colData(ref_)[anno_label]))
if(!col_file=="nan"){col_in=readRDS(col_file)}
if(col_file=="nan"){
  df_cols=colData(ref_)[c(anno_label,"cols")] %>% as.data.frame() %>% distinct()
  names(df_cols)[1]="anno_label"
  col_in=df_cols$cols
  names(col_in)=df_cols$anno_label
}

#read seurat obj  ------------------------------------------------------------------------------------- ----
cat("Reading obj ...\n")

obj_QC=readRDS(rds_in)
print(obj_QC)

# #subset object for test  ------------------------------------------------------------------------------------- ----
if(do_subset){
  set.seed(10)
  df_barcode=get_barcode_df(obj_QC) %>% sample_n(n_max)
  obj_QC=subset(obj_QC,cells=df_barcode$barcode)
}

#redo cluters  ------------------------------------------------------------------------------------- ----
if(do_clustering){
  cat("Running clustering ...\n")
  obj_QC=FindNeighbors(obj_QC)
  obj_QC=FindClusters(obj_QC,resolution = 2)
}

DimPlot(obj_QC,group.by = "seurat_clusters")
ggsave(paste0(dir_out,"dimPlot.seurat_clusters",".png"),width=12,height = 10)

#annotation  ------------------------------------------------------------------------------------- ----
cat("Running SingleR ...\n")

singleR_out_=run_singleR(inObj = obj_QC,ref_,anno_label,transform = transform)

obj_QC=dimPlot_fromDF(objIn = obj_QC,df = singleR_out_$df_singleR_out,
                            var = "labels",
                            label = paste0(sample,".",label,".Celllevel"),
                            cols_in = col_in,
                            dir_out = dir_out)


obj_QC=dimPlot_fromDF(objIn = obj_QC,df = singleR_out_$df_singleR_out,
                            var = "labels_cluster",
                            label = paste0(sample,".",label,".Clusterlevel"),
                            cols_in = col_in,
                            dir_out = dir_out)

file_singleR_out=paste0(dir_out,sample,".",label,".SingleR_out.rds")

cat(paste0("Output singleR file: ",file_singleR_out,"\n"))

saveRDS(singleR_out_,paste0(file_singleR_out))

#get stringent labels  ------------------------------------------------------------------------------------- ----
# df_label_stringent=get_barcode_df(obj_QC) %>% 
#   left_join(get_embedding_df(obj_QC,"umap")) %>% 
#   left_join(get_stringent_label(singleR_out_$out_singleR,0.75))
# 
# write_tsv(df_label_stringent,paste0(dir_out,"df_label_stringent",".tsv"))
# 
# ggplot(df_label_stringent[!is.na(df_label_stringent$label_stringent),],aes(x=UMAP_1,y=UMAP_2)) + 
#   geom_point(aes(color=label_stringent)) + theme_bw() + ggtitle("_Stringent") +
#   scale_color_manual(values = col_in[names(col_in) %in% df_label_stringent$label_stringent]) 
# 
# ggsave(paste0(dir_out,"SingleR_annotation_","_Stringent",".png"),width=10,height = 10)
# ggsave(paste0(dir_out,"SingleR_annotation_","_Stringent",".pdf"),width=10,height = 10)

#plot delta  ------------------------------------------------------------------------------------- ----
cat("Running SingleR QC ...\n")

SingleR::plotDeltaDistribution(singleR_out_$out_singleR)
ggsave(paste0(dir_out,sample,".",label,".DeltaDistribution_cellLevel",".png"),width=15,height = 12)

SingleR::plotDeltaDistribution(singleR_out_$out_singleR_cluster)
ggsave(paste0(dir_out,sample,".",label,".DeltaDistribution_ClusterLevel",".png"),width=15,height = 12)

#output  ------------------------------------------------------------------------------------- ----

# saveRDS(obj_QC,file=paste0(dir_out,sample,".",label,".Seurat.rds"))

cat("Done SingleR\n")
























