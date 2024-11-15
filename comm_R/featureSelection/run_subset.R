# gc(rm(list=ls()))
rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

library (SummarizedExperiment, quietly = T)
library(dplyr)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# 
# obj_rds=              "../tsne/2042_junctions/obj_2042_junctions.rds"
# dataset_label=        "2042_junctions"
# subset_label=         "1st"
# feature_subset_file=  "2042_junctions/Boruta_summarise/Boruta_2042_junctions_ConfirmedTentative.list"
# dir_out=              "2042_junctions/subset/"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=              args[1]
dataset_label=        args[2]
feature_subset_file=  args[3]
subset_label=         args[4]
dir_out=              args[5]


features=NULL
#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
print("Reading Obj")
obj_=readRDS(obj_rds)

#get features ----
print("Get feature names")

features_in=readLines(feature_subset_file)

df_feature=data.frame(features=features_in,stringsAsFactors = F) %>% 
    arrange(features) %>% 
    mutate(obs=1:n())

print(paste0("Totall included number of features: ",nrow(df_feature)))

# row.names(obj_$SE)

#create output files ---------------
  obj_x=obj_
  obj_x$SE=obj_x$SE[df_feature$features,]
  saveRDS(obj_x,paste0(dir_out,"obj_",dataset_label,"_subset",subset_label,".rds"))


