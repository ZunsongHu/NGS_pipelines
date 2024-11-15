# rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

library(SummarizedExperiment, quietly = T)
library(caret,quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# 
# rds_obj="2707/sampleSelection/5/get_topMAD/obj.2707.topMAD5000.rds"
# file_umap="2655/sampleSelection/1/umap/TopMAD5000/2655.TopMAD5000.FeatureN1200.NeighborsN10.umap.tsv"
# dir_out="test/"
# dataset_label="2655"
# feature_panel="TopMAD5000"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_obj=      args[1]
file_umap=    args[2]
dir_out=      args[3]
dataset_label=args[4]
feature_panel=args[5]

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read data ----
obj_=readRDS(rds_obj)

df_umap=read_tsv(file_umap)

df_merge=get_features_df(obj_in = obj_,features = c("COH_sample","diag"))  %>% left_join(df_umap)

NeighborN=unique(df_merge$n_neighbors)
FeatureN=unique(df_merge$FeatureN)
#run KNN ----
names(df_merge)
df_merge$KNN_pred=knn_pred_one(indata = df_merge,var_y = "diag",var_x = c("uMAP_1","uMAP_2"),KNN_K = 3)

file_output=paste0(dir_out,dataset_label,".",feature_panel,".FeatureN",FeatureN,".NeighborsN",NeighborN,".umap_KNN.tsv")

cat(paste0(file_output,"\n"))
write_tsv(df_merge,file_output)















































