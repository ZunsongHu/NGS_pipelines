#set ---------------------------
rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
source("/home/zuhu/bin/R/plots/F_plots.R")


#testing parameters ----
setwd("/scratch/zuhu/project/BALL/feature_selection/")

file_list="2391_totalRNA/get_errorRate/df_error.list"
datasetLabel="2391_totalRNA"
dir_out="test/"

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

#read dfs for error rate ----

file_list=readLines(file_list)

df_error=bind_rows(lapply(file_list, function(file_one){
  read_tsv(file_one)
}))
  
write_tsv(df_error,paste0(dir_out,"df_error.",datasetLabel,".FeaturePanel.",feature_panel,".tsv"))

#get error plot ----

gg_linePlot(df_error,x = "featureN",y = "error_rate",var_col = "diag",cols = subtypeCol,split.by="diag")
ggsave(paste0(dir_out,"Plot_error.",datasetLabel,".FeaturePanel.",feature_panel,".pdf"),width=40,height=20)






























  





