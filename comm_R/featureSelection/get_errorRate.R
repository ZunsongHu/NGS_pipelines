#set ---------------------------
rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# 
# file_rds_obj="../tsne/2391/obj_2391_totalRNA.rds"
# datasetLabel="2391_totalRNA"
# feature_panel="totalRNA"
# featureN=as.numeric("60")
# dir_out="test/"
  
#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_rds_obj=       args[1]
datasetLabel=       args[2]
feature_panel=      args[3]
featureN=           as.numeric(args[4])
dir_out=            args[5]

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read object ------
cat("Reading obj...")
obj_=readRDS(file_rds_obj)

#run prediction ----
obj_=run_PhenoGraph(obj_in = obj_,feature_panel = feature_panel,variable_n = featureN,neighbor_k = 3)
df_pred=obj_$PhenographPred[c("COH_sample","diag","diag_pred")] 

#get error rate ----
df_error=df_pred %>% 
  group_by(diag) %>% 
  mutate(N_EachSubtype=n(),
         success_pred=ifelse(diag==diag_pred,"Yes","No")) %>% 
  group_by(diag,success_pred) %>% 
  mutate(N_successPredEachSubtype=n()) %>% 
  ungroup() %>% 
  filter(success_pred=="Yes") %>% 
  select(diag,N_EachSubtype,N_successPredEachSubtype) %>% 
  distinct() %>% 
  mutate(N_failPredEachSubtype=N_EachSubtype-N_successPredEachSubtype,
         error_rate=round(N_failPredEachSubtype/N_EachSubtype,4),
         datasetLabel=datasetLabel,
         feature_panel=feature_panel,
         featureN=featureN)

write_tsv(df_error,paste0(dir_out,"df_error.",datasetLabel,".",feature_panel,'.featureN',featureN,".tsv"))

































