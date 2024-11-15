# rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
# source("/home/zuhu/bin/R/plots/F_plots.R")

library(SummarizedExperiment, quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# 
# obj_rds=            "../tsne/2042/obj_2042.rds"
# label=              "2042"
# boruta_confirmed=   "2042/Boruta_summarise/Boruta_2042_ConfirmedList.tsv"
# dir_out=            "test/"


#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=            args[1]
label=              args[2]
boruta_confirmed=   args[3]
dir_out=            args[4]

confirmed_times_cutoff=3

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
cat("Read obj...\n")
obj_=readRDS(obj_rds)

#read boruta results ----
df_boruta_confirmed=read_tsv(boruta_confirmed) %>% arrange(rank_median)

df_boruta_confirmed1=df_boruta_confirmed %>% 
  filter(confirmed_times >= confirmed_times_cutoff)  %>% 
  mutate(source=label)

df_boruta_confirmed1_featured=df_boruta_confirmed1 %>% 
  filter(!mean3SD_picked=="")

#remove corr features ----
obj_[[paste0("confirmed")]]=df_boruta_confirmed1$feature_id
obj_[[paste0("confirmed_noCorr")]]=get_nonCorr_genes(obj_in = obj_,feature_panel = paste0("confirmed"),featureRank = df_boruta_confirmed1$rank_median,highCorrCutoff = 0.9)

obj_[[paste0("confirmed_featured")]]=df_boruta_confirmed1_featured$feature_id
obj_[[paste0("confirmed_featured_noCorr")]]=get_nonCorr_genes(obj_in = obj_,feature_panel = paste0("confirmed_featured"),featureRank = df_boruta_confirmed1_featured$rank_median,highCorrCutoff = 0.75)

#Get top MAD features ----
obj_mad=obj_

obj_mad$MAD_tmp=get_variable_genes(obj_in = obj_mad,N_genes = 2000)

obj_mad$MAD_tmp_noCorr=get_nonCorr_genes(obj_in = obj_mad,feature_panel ="MAD_tmp",featureRank = 1:length(obj_$MAD_tmp),highCorrCutoff = 0.75)

obj_mad$SE=obj_mad$SE[obj_mad$MAD_tmp_noCorr,]

obj_$variable_genes=get_variable_genes(obj_in = obj_mad,N_genes = 1500)

#Svae output -----
obj_out=obj_
# obj_out$SE=NULL
obj_out$df_boruta_confirmed=df_boruta_confirmed1
obj_out$df_boruta_confirmed_featured=df_boruta_confirmed1_featured

saveRDS(obj_out,paste0(dir_out,"obj_borutaConfirmedFeatures_",label,".rds"))









