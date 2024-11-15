#rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
library (SummarizedExperiment, quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# 
# obj_rds="2655/sampleSelection/1/get_topMAD/obj.2655.topMAD5000.rds"
# dir_out="test/"
# dataset_label="2655"
# feature_panel="TopMAD5000"
# variable_n=800

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=          args[1]
dir_out=          args[2]
dataset_label=    args[3]
feature_panel=    args[4]
variable_n=       as.numeric(args[5])

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
cat("Reading obj...\n")
obj_=readRDS(obj_rds)

#run tsne ----
obj_PCA=run_PCA(obj_in = obj_,dims=100,variable_n=variable_n,feature_panel=feature_panel,out_label="pca")

df_PCA=obj_PCA$pca %>% mutate(FeatureN=variable_n)
df_PCA$COH_sample=row.names(df_PCA)

cat(paste0(dir_out,dataset_label,".",feature_panel,".FeatureN",variable_n,".PCA.tsv\n"))
write_tsv(df_PCA,paste0(dir_out,dataset_label,".",feature_panel,".FeatureN",variable_n,".PCA.tsv"))
















































