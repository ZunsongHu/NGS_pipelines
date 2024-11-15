# rm(list=ls())

library(stringr)
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

# #temp parameters ----
setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq/ScPCA_Portal/")

file_obj="out_raw/SCPCL000603/Seurat/SCPCL000603.SeuratObj.QC.rds"
feature_panel_label="CRLF2"
file_feature="CRLF2"
dir_out="test/"
file_singleR_celltype="out_raw/SCPCL000603/BALL_subtyping/obj_singleR_celltype.SCPCL000603.rds"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_obj=args[1]
feature_panel_label=args[2]
file_feature=args[3]
dir_out=args[4]
file_singleR_celltype=args[5]

id=gsub(".SeuratObj.QC.rds","",basename(file_obj))

#create dir out ----
dir_out_=paste0(dir_out,feature_panel_label,"/")
if(!dir.exists(dir_out)){dir.create(dir_out)}
if(!dir.exists(dir_out_)){dir.create(dir_out_)}

#read obj ----
obj_=readRDS(file_obj)

#get key genes matrix ----
feature_genes=readLines(file_feature)

df_matrix=as.data.frame(as.matrix(t(GetAssayData(obj_))))
df_matrix$barcode=row.names(df_matrix)
genes_intersect=colnames(df_matrix)[toupper(colnames(df_matrix)) %in% feature_genes]
df_matrix_gene=df_matrix[c("barcode",genes_intersect)]

feature_genes_in=genes_intersect

genes_un=feature_genes[!feature_genes %in% genes_intersect]
if(length(genes_un)>0){
  df_matrix_count=as.data.frame(as.matrix(t(GetAssayData(obj_,assay = "RNA",slot = "count"))))
  df_matrix_count$barcode=row.names(df_matrix_count)
  genesUn_intersect=genes_un[genes_un %in% toupper(colnames(df_matrix_count))]
  
  if(length(genesUn_intersect) >=1){
    df_matrix_gene=df_matrix_gene %>% left_join(df_matrix_count[c("barcode",genesUn_intersect)])
  }
  feature_genes_in=c(genes_intersect,genesUn_intersect)
  
}

#get umap ----
df_umap=get_barcode_df(obj_) %>% 
  left_join(get_embedding_df(obj_,"umap")) %>%
  left_join(df_matrix_gene) 

#celltype ---------------
if(!(file_singleR_celltype=="NULL" | is.null(file_singleR_celltype))){
  obj_singleR=readRDS(file_singleR_celltype)
  
  df_singleR=get_barcode_df(obj_) %>%
    left_join(obj_singleR$df_singleR_out[c("barcode","labels")])

  obj_[["SingleR_celltype"]]=df_singleR$labels

  names(df_singleR)[2]="SingleR_celltype"
  df_celltype=df_umap %>% left_join(df_singleR)
}

#cell cycle ----
df_cellcycle=df_celltype %>% left_join(get_metadataVar(inObj = obj_,var = c("Phase")))

#save output file -----------------
write_tsv(df_cellcycle,paste0(dir_out_,"df_features.",id,".",feature_panel_label,".tsv"))

#Draw dot plot ----
save_dot_plot=function(obj=obj_,group="seurat_clusters"){
  DotPlot(object = obj, features = feature_genes_in,group.by = group)
  ggsave(paste0(dir_out_,"dotPlot.",id,".",group,".png"),width=20,height = 10)
  ggsave(paste0(dir_out_,"dotPlot.",id,".",group,".pdf"),width=20,height = 10)
}

save_heatmap=function(obj=obj_,group="seurat_clusters"){
  df_sct=GetAssayData(obj,slot = "scale",assay = "SCT")
  if(any(feature_genes_in %in% row.names(df_sct))){
    DoHeatmap(object = obj, features = feature_genes_in,group.by = group)
    ggsave(paste0(dir_out_,"heatmap.",id,".",group,".png"),width=20,height = 10)
    ggsave(paste0(dir_out_,"heatmap.",id,".",group,".pdf"),width=20,height = 10)
  } else {
    message("No feature expressed high enough to do Heatmap!")
  }
}

save_VlnPlot=function(obj=obj_,group="seurat_clusters"){
  VlnPlot(object = obj, features = feature_genes_in,group.by = group)
  ggsave(paste0(dir_out_,"VlnPlot.",id,".",group,".png"),width=20,height = 10)
  ggsave(paste0(dir_out_,"VlnPlot.",id,".",group,".pdf"),width=20,height = 10)
}

save_dot_plot(obj=obj_,group="SingleR_celltype")
save_dot_plot(obj=obj_,group="Phase")

save_heatmap(obj=obj_,group="SingleR_celltype")
save_heatmap(obj=obj_,group="Phase")

save_VlnPlot(obj=obj_,group="SingleR_celltype")
save_VlnPlot(obj=obj_,group="Phase")

#Feature bar plot ----
save_barPlot=function(df_in=df_cellcycle,var=feature_genes_in[1],group="Phase"){
  
  df_in["var"]=unlist(df_in[var])
  df_in["group"]=unlist(df_in[group])
  
  df_in=df_in %>% mutate(var_g=ifelse(var>=median(var),"1High","0Low"))
  
  if(length(levels(as.factor(df_in$var_g)))<2){
    message(paste0("Group based on median of ",feature, " failed! No barPlot for ",feature,"!"))
  }
  
  if(length(levels(as.factor(df_in$var_g)))>=2){
    df_in1=df_in %>% select(var_g,group) %>% group_by(var_g) %>% mutate(n_g=n()) %>% group_by(var_g,group) %>% mutate(n=n(),per=n/n_g) %>% distinct() %>% arrange(var_g,group)
    
    gg_barplot(data_in = df_in1,x = "var_g",y = "per",group.by = "group",position = "stack",x_lab = paste0(var,"(group by median)"),y_lab="Percentages") +
      theme(legend.position = "right") +
      guides(fill = guide_legend(title.position = "top",title.hjust = 0.5,title = group)) +
      labs(title = paste0(var,"(group by median) vs ",group))
    
    ggsave(paste0(dir_out_,"BarPlot.",id,".",var,"_",group,".png"),width=6,height = 5)
    ggsave(paste0(dir_out_,"BarPlot.",id,".",var,"_",group,".pdf"),width=6,height = 5)
    
    names(df_in1)[1:2]=c(var,group)
    write_tsv(df_in1,paste0(dir_out_,"df_percentage_",id,".",var,"_",group,".tsv"))
  }
}

for (feature in feature_genes_in){
save_barPlot(df_in=df_cellcycle,var=feature,group="Phase" )
}











