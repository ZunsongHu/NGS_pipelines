# rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

library(SummarizedExperiment, quietly = T)

#testing parameters ----
setwd("/scratch/zuhu/project/BALL/feature_selection/")

smote=            "2042_middleGroup_mRNA/smote/smote_2042_middleGroup_mRNA_diag_EachGroupN10.rds"
label=            "2042_middleGroup_mRNA"
boruta_out=       "2042_middleGroup_mRNA/Boruta_summarise/Boruta_2042_middleGroup_mRNA_rawOut.tsv"
smote_N=          10
dir_out=          "test/"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

smote=        args[1]
label=        args[2]
boruta_out=   args[3]
smote_N=      as.numeric(args[4])
dir_out=      args[5]

#create dir out ----
dir_out_umap=paste0(dir_out,"umap/")
dir_out_tsne=paste0(dir_out,"tsne/")

if(!dir.exists(dir_out)){dir.create(dir_out)}
if(!dir.exists(dir_out_umap)){dir.create(dir_out_umap)}
if(!dir.exists(dir_out_tsne)){dir.create(dir_out_tsne)}

#read boruta ----
df_boruta=read_tsv(boruta_out) %>% 
  filter(EachGroupN==smote_N & decision=="Confirmed") %>% 
  arrange(desc(meanImp))

#read smote ----
cat("Read smote...\n")
smote_=readRDS(smote)

df_info=smote_["diag"]

matrix_=smote_ %>% mutate(diag=NULL) %>% 
  t() %>% as.matrix()

#get obj ----
obj_=list(SE=SummarizedExperiment(assays=list(vst=matrix_), colData = df_info))

obj_[["confirmed"]]=df_boruta$feature_id

#run umap ----
variable_n_list=c(55,100,200,500,800,1000,1500,2000,2500,length(df_boruta$feature_id))
variable_n_list=variable_n_list[variable_n_list<=length(df_boruta$feature_id)]

# variable_n_list=55
# variable_n=55

for(variable_n in variable_n_list){
  obj_=run_umap(obj_in = obj_,n_neighbors = 10,dims = 2,variable_n = variable_n,out_label = "umap",feature_panel = "confirmed")
  
  write_tsv(obj_$umap,paste0(dir_out_umap,"umap.smoteN",smote_N,".FeatureN",variable_n,".tsv"))
  
  draw_DimPlot(obj_in = obj_,group.by = 'diag',reduction = "umap")
  ggsave(paste0(dir_out_umap,"umap.smoteN",smote_N,".FeatureN",variable_n,".pdf"),width = 11,height = 9)
  ggsave(paste0(dir_out_umap,"umap.smoteN",smote_N,".FeatureN",variable_n,".png"),width = 11,height = 9)
}

#run tsne ----

for(variable_n in variable_n_list){
  
  perplexity_one=ifelse(smote_N>=20,20,smote_N)
  
  obj_=run_tsne(obj_in = obj_,perplexity_one = perplexity_one,dims = 2,variable_n = variable_n,out_label = 'tsne',feature_panel = "confirmed",initial_dims = 49)
  
  write_tsv(obj_$tsne,paste0(dir_out_tsne,"tsne.smoteN",smote_N,".FeatureN",variable_n,".tsv"))
  
  draw_DimPlot(obj_in = obj_,group.by = 'diag',reduction = "tsne")
  print(paste0(dir_out_tsne,"tsne.smoteN",smote_N,".FeatureN",variable_n,".pdf"))
  ggsave(paste0(dir_out_tsne,"tsne.smoteN",smote_N,".FeatureN",variable_n,".pdf"),width = 11,height = 9)
  ggsave(paste0(dir_out_tsne,"tsne.smoteN",smote_N,".FeatureN",variable_n,".png"),width = 11,height = 9)
}










