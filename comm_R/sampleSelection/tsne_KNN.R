
gc(rm(list=ls()))

source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

library(SummarizedExperiment, quietly = T)
library(caret,quietly = T)
#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_tsne=    args[1]
dataset_label=args[2]
file_diag=    args[3]
dir_out=      args[4]


#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# file_tsne="sample_selection/1/tsne/1.perplexityN10.FeatureN500.tsv"
# dataset_label="1"
# file_diag="id/df_info_2955_1.tsv"
# dir_out="test/"

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read data ----
df_diag=read_tsv(file_diag)

df_tsne=read_tsv(file_tsne)

df_merge=df_tsne %>% select(COH_sample,tSNE_1,tSNE_2,"perplexityN","FeatureN" ) %>% left_join(df_diag %>% select(COH_sample,diag))

perplexityN=unique(df_merge$perplexityN)
FeatureN=unique(df_merge$FeatureN)
#run KNN ----
df_merge$KNN_pred=knn_pred_one(indata = df_merge,var_y = "diag",var_x = c("tSNE_1","tSNE_2"),KNN_K = 3)

write_tsv(df_merge,paste0(dir_out,dataset_label,".KNN.perplexityN",perplexityN,".FeatureN",FeatureN,".tsv"))











































