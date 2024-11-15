
rm(list=ls())

source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

library(SummarizedExperiment, quietly = T)
library(caret,quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# 
# file_pca="2655/sampleSelection/1/PCA/TopMAD5000/2655.TopMAD5000.FeatureN800.PCA.tsv"
# rds_obj="2655/sampleSelection/1/get_topMAD/obj.2655.topMAD5000.rds"
# dataset_label="2655"
# feature_panel="TopMAD5000"
# used_PC_n=50
# dir_out="test/"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_obj=      args[1]
file_pca=     args[2]
dir_out=      args[3]
dataset_label=args[4]
feature_panel=args[5]
used_PC_n=    as.numeric(args[6])


#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read data ----
obj_=readRDS(rds_obj)

df_pca=read_tsv(file_pca)

df_merge=get_features_df(obj_in = obj_,features = c("COH_sample","diag"))  %>% left_join(df_pca)

#run KNN ----
FeatureN=unique(df_merge$FeatureN)

var_x_in=paste0("PC",1:used_PC_n)

df_merge$KNN_pred=knn_pred_one(indata = df_merge,var_y = "diag",var_x = var_x_in,KNN_K = 5)

df_out=df_merge  %>% select(COH_sample,diag,KNN_pred,FeatureN) %>% mutate(PC_N=used_PC_n)


file_output=paste0(dir_out,dataset_label,".",feature_panel,".FeatureN",FeatureN,".PCN",used_PC_n,".PCA_KNN.tsv")

cat(paste0(file_output,"\n"))
write_tsv(df_out,file_output)











































