# options(future.globals.maxSize = 10*1024 * 1024^2) #10GB
options(future.globals.maxSize = 150*1024 * 1024^2) #150GB

rm(list=ls())
library(ggplot2)
library(Seurat)

source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

#temp parameters
# setwd("/scratch/zuhu/tmp/SC_CMV/")
# 
# rds_in=     "out_raw/integration/threesample/merge/obj_merge_1.rds"
# rds_out=     "out_raw/integration/threesample/merge/obj_merge_1_batchCorrected.rds"

#get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=args[1]
rds_out=args[2]

#create dir_out ----
dir_out=dirname(rds_out)
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
obj_=readRDS(rds_in)

obj_corrected=batch_correction_seurat(inObj = obj_,var_batch="group",n_varaible=3000,method="CCA")

# obj_corrected=FindVariableFeatures(obj_corrected,nfeatures = 2000)

obj_corrected=RunPCA(obj_corrected)

obj_corrected=RunUMAP(obj_corrected,dims = 1:50)

saveRDS(obj_corrected,rds_out)

