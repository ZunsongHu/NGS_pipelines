# gc(rm(list=ls()))
rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

library (SummarizedExperiment, quietly = T)
library(dplyr)

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=          args[1]
dataset_label=    args[2]
feature_N_eachFile=  args[3]
pieces_N=          as.numeric(args[4])
dir_out=          args[5]

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# 
# obj_rds=          "../tsne/2042_junctions/obj_2042_junctions.rds"
# dataset_label=    "2042_junctions"
# feature_N_eachFile="NULL"
# pieces_N=         "12"
# dir_out=          "2042_junctions/split/"

random_split=T

features=NULL
#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
print("Reading Obj")
obj_=readRDS(obj_rds)

#get features ----
print("Get feature names")
if(is.null(features)){features_in=rownames(obj_$SE)}
if(!is.null(features)){features_in=features}
# features_in=features_in[1:100]

set.seed(10)
if(random_split){
  df_feature=data.frame(features=features_in,stringsAsFactors = F) %>% 
    mutate(random_number=runif(n())) %>% 
    arrange(random_number) %>% 
    mutate(obs=1:n(),random_number=NULL)
}

if(!random_split){
  df_feature=data.frame(features=features_in,stringsAsFactors = F) %>% 
    arrange(features) %>% 
    mutate(obs=1:n())
}


#get number of pieces and number of feature in each piece
if(!feature_N_eachFile=="NULL" & !pieces_N=="NULL"){stop("Only one parameter (feature_N_eachFile or pieces_N) is allowed")}

if(feature_N_eachFile=="NULL" & !pieces_N=="NULL"){
  pieces_N=as.numeric(pieces_N)
  feature_N_eachFile=ceiling(nrow(df_feature)/pieces_N)
}

if(!feature_N_eachFile=="NULL" & pieces_N=="NULL"){
  feature_N_eachFile=as.numeric(feature_N_eachFile)
  pieces_N=ceiling(nrow(df_feature)/feature_N_eachFile)
}

cat(paste0("Total pieces: ",pieces_N,"; Number of features in each piece: ",feature_N_eachFile,"\n"))

#get output index  ---------------
df_feature$pieces_i=NA
for (i in 1:pieces_N){
  df_feature$pieces_i[df_feature$obs > feature_N_eachFile*(i-1) & df_feature$obs <= feature_N_eachFile*i]=i
}

# df_feature %>% group_by(pieces_i) %>% mutate(max=max(obs)) %>% select(max) %>% distinct()

#create output files ---------------
for (i in 1:pieces_N){
  cat(paste0("Creating files: ",i,"\n"))
  obj_x=obj_
  obj_x$SE=obj_x$SE[df_feature$features[df_feature$pieces_i==i],]
  saveRDS(obj_x,paste0(dir_out,"obj_",dataset_label,"_split",i,".rds"))
}


