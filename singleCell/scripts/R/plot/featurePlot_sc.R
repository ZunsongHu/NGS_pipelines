# rm(list=ls())

library(stringr)
# source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")
# source("/home/zuhu/bin/R/plots/F_plots.R")

# functions ------------------------------------------------------------------------------------------- ----

draw_feature_plot_fromDF=function(objIn,df,varname,label,cols_in,dir_out){
  
  objIn[[label]]=unlist(df[varname])
  
  DimPlot(objIn,group.by = label,label=T,repel = T,
          cols = cols_in[names(cols_in) %in% unlist(objIn@meta.data[label])])
  ggsave(paste0(dir_out,label,".png"),width=16,height = 9)
  ggsave(paste0(dir_out,label,".pdf"),width=16,height = 9)
  
  objIn
}


draw_umap_continous=function(inObj,var,outfile_prefix="2.data_temp/test/VlnPlot_qc"){
  dir_name=dirname(dirname(outfile_prefix));if(!dir.exists(dir_name)){dir.create(dir_name)}
  dir_name=dirname(outfile_prefix);if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  FeaturePlot(inObj, features = var)
  ggsave(filename = paste0(outfile_prefix,".png"),width = 13,height = 9)
  ggsave(filename = paste0(outfile_prefix,".pdf"),width = 13,height = 9)
}



#temp parameters ----
setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq/GuLab/")

file_obj="out_sc/COH005637D1/Seurat/COH005637D1.SeuratObj.QC.rds"
feature_panel_label="CRLF2"
file_feature="CRLF2"
dir_out="test/"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_obj=args[1]
feature_panel_label=args[2]
file_feature=args[3]
dir_out=args[4]

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


#seurat_clusters and gene matrix ----
df_umap=get_barcode_df(obj_) %>% left_join(get_embedding_df(obj_,"umap")) 

df_out=df_umap %>% left_join(df_matrix_gene)

#save EXP matrix -----------------
write_tsv(df_out,paste0(dir_out_,"df_features_",feature_panel_label,".tsv"))

#draw scatter plot ----
# df_in=df_out
# x="UMAP_1"
# y="UMAP_2"
# var_col="AIM2"

for(feature_one in feature_genes_in){
  gg_featurePlot(df_in = df_out,x="UMAP_1",y="UMAP_2",var_col = feature_one,size = 0.3)
  ggsave(paste0(dir_out_,"featurePlot.",id,".",feature_one,".png"),width=7,height = 7)
  ggsave(paste0(dir_out_,"featurePlot.",id,".",feature_one,".pdf"),width=7,height = 7)
}






