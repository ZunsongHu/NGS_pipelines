# usage
# Rscript_seurat5 /net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/singleCell/scripts/R/plot/dimPlot_sc.R \
# analysis/obj/obj_seurat.harmony.10samples.rds \
# sample \
# analysis/dimPlot/group/dimPlot.sample.pdf \
# F F \
# 8 6



# options(future.globals.maxSize = 150*1024 * 1024^2) #50GB
# options(scipen = 999)

# options(warn=-1)
# suppressPackageStartupMessages({
#   library(Seurat)
#   library(dplyr)
#   library(stringr)
#   library(ggplot2)
#   library(scCustomize)
# })
# options(warn=0)


# functions  ------------------------------------------------------------------------------------- ----
get_barcode_df=function(object_in){
  data.frame(barcode=row.names(object_in@meta.data),stringsAsFactors = F)
}

# temp parameters ------------------------------------------------------------------------------------- ----
# setwd("/net/nfs-irwrsrchnas01/labs/zgu_grp/Group/Grp_AncaPasca/omics/")
# 
# file_rds="analysis/obj/obj_seurat.harmony.5samples.rds"
# var="sample"
# pdf_out="analysis/dimPlot/sample/dimPlot.sample.pdf"
# 
# label_figure='F';repel_figure='F'
# w_umap=7;h_umap=6

# parameters ------------------------------------------------------------------------------------- ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_rds=args[1]
var=args[2]
pdf_out=args[3]

label_figure=args[4];repel_figure=args[5]
w_umap=args[6];h_umap=args[7]


# other parameters ------------------------------------------------------------------------------------- ----
shuffle=T

# pharse parameters ------------------------------------------------------------------------------------- ----

dir_out=paste0(dirname(pdf_out),"/")
if(!dir.exists(dir_out)){dir.create(dir_out,recursive = T)}
prefix=gsub(".pdf","",pdf_out)

label_figure=ifelse(label_figure=="F",F,T)
repel_figure=ifelse(repel_figure=="F",F,T)

w_umap=as.numeric(w_umap);h_umap=as.numeric(h_umap)

# read obj  ------------------------------------------------------------------------------------- ----
cat("reading obj...\n")

obj_=readRDS(file_rds)


# get umap   ------------------------------------------------------------------------------------- ----
DimPlot(obj_,group.by = var,label = label_figure,repel = repel_figure,shuffle = shuffle)
ggsave(pdf_out,width = w_umap,height = h_umap)
ggsave(paste0(prefix,".png"),width = w_umap,height = h_umap)




cat("Done umap!\n")




# functions ------------------------------------------------------------------------------------------- ----
draw_DimPlot=function(inObj,df_feature,var_key,label,file_out_prefix){
  if (!all(row.names(inObj@meta.data)==unlist(df_feature[1]))){stop("barcode not equal")}
  
  inObj[[label]]=unlist(df_feature[var_key])
  DimPlot(inObj,group.by = label,label=T,repel = T,cols=DiscretePalette(length(unique(unlist(df_feature[var_key])))))
  
  ggsave(paste0(file_out_prefix,".png"),width=10,height = 10)
  ggsave(paste0(file_out_prefix,".pdf"),width=10,height = 10)
  inObj
}

draw_umap_mini=function(obj,var_group,label,dir_out,width=15,height=11){
  
  n_levels=length(unique(unlist(get_metadataVar(obj,var_group)[,2])))
  
  DimPlot(obj,group.by = var_group,label=T,repel = T,cols=DiscretePalette(n_levels))
  ggsave(paste0(dir_out,"umap_",label,"_",var_group,".pdf"),width=width,height = height)
  ggsave(paste0(dir_out,"umap_",label,"_",var_group,".png"),width=width,height = height)
  
}


draw_umap_class=function(inObj,var,outfile_prefix="2.data_temp/test/VlnPlot_qc"){
  dir_name=dirname(dirname(outfile_prefix));if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  dir_name=dirname(outfile_prefix);if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  var_value=unlist(inObj@meta.data[var])
  
  Idents(inObj)=var_value
  
  col_in=unique(subtypeCol)[2:length(unique(subtypeCol))][1:length(unique(var_value))]
  
  DimPlot(inObj,cols = col_in)
  ggsave(filename = paste0(outfile_prefix,".png"),width = 13,height = 9)
  ggsave(filename = paste0(outfile_prefix,".pdf"),width = 13,height = 9)
}



dimPlot_fromDF=function(objIn = obj_QC,df = singleR_out_$df_singleR_out,
                        var = "labels",
                        label = paste0(sample,".",label,".Celllevel"),
                        cols_in = col_in,
                        dir_out = dir_out,w=12,h=9,
                        label_figure=F,repel_figure=F
){
  
  if(!all(colnames(obj_QC)==df$barcode)){stop("barcodes of obejct and df for var are not matched")}
  
  objIn[[label]]= unname(unlist(df[var]))
  
  DimPlot(objIn,group.by = label,label=label_figure,repel = repel_figure,
          cols = cols_in[names(cols_in) %in% unlist(objIn@meta.data[label])])
  ggsave(paste0(dir_out,label,".png"),width=w,height = h)
  ggsave(paste0(dir_out,label,".pdf"),width=w,height = h)
  objIn
}



























