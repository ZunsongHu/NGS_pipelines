gc(rm(list=ls()))
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

# library(SummarizedExperiment, quietly = T)
# library(caret,quietly = T)
#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

vst_matrix=   args[1]
file_diag=    args[2]
label=        args[3]
dir_out=      args[4]

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/2295/")
# vst_matrix=   "/scratch/zuhu/project/BALL/tsne/2863/vst_2863.rds"
# file_diag=      "tsneKNN/2295.tsneKNN_Matched.tsv"
# label=      "test"
# dir_out=      "test/"

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

vst_=readRDS(vst_matrix)
row.names(vst_)=vst_$feature

df_diag=read_tsv(file_diag)
vst_out=vst_[c("feature",df_diag$COH_sample)]

row.names(vst_out)=vst_out$feature

saveRDS(vst_out,paste0(dir_out,label,".vst.rds"))

















