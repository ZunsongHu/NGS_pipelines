# rm(list=ls())

source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
library (SummarizedExperiment, quietly = T)
library(Rphenograph)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# 
# obj_rds="2707/sampleSelection/1/get_topMAD/obj.2707.topMAD5000.rds"
# dir_out="test/"
# dataset_label="2707"
# feature_panel="TopMAD5000"
# neighbor_k=10
# variable_n=800

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=      args[1]
dir_out=      args[2]
dataset_label=args[3]
feature_panel=args[4]
neighbor_k=   as.numeric(args[5])
variable_n=   as.numeric(args[6])

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
cat("Reading obj...\n")
obj_=readRDS(obj_rds)

#run tsne ----
obj_out=run_PhenoGraph(obj_in = obj_,variable_n = variable_n,neighbor_k = neighbor_k,feature_panel = feature_panel)
df_PhenographPred=obj_out$PhenographPred

write_tsv(df_PhenographPred,paste0(dir_out,dataset_label,".",feature_panel,".FeatureN",variable_n,".NeighborN",neighbor_k,".PhenographPred.tsv"))
























