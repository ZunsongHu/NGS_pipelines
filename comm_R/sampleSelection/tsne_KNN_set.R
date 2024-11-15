gc(rm(list=ls()))
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

# library(SummarizedExperiment, quietly = T)
# library(caret,quietly = T)
#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_KNN_list=    args[1]
dataset_label=args[2]
dir_out=      args[3]
cutoff_diag_n= as.numeric(args[4])


#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# file_KNN_list=    "2280/tsneKNN/2280.tsneKNN.list"
# dataset_label="test"
# dir_out=      "test"

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read files list ----
files=readLines(file_KNN_list)

#get kept level ----
df_one=read_tsv(files[1])
diag=df_one$diag
diag_level_keep=(as.data.frame(table(diag)) %>% mutate(diag=as.character(diag)) %>% arrange(Freq) %>% filter(Freq <= cutoff_diag_n))$diag

#get output ----
df_out=bind_rows(lapply(files, function(file_one){
  read_tsv(file_one)
}))

df_filter_list=df_out %>% filter(!diag == KNN_pred) %>% select(COH_sample) %>% distinct()

df_out1=df_out %>% filter(
  (!COH_sample %in% df_filter_list$COH_sample) | diag %in% diag_level_keep
) %>% select(COH_sample,diag) %>% distinct()

write_tsv(df_out,paste0(dir_out,dataset_label,".tsneKNN_set.tsv"))
write_tsv(df_out1,paste0(dir_out,dataset_label,".tsneKNN_Matched.tsv"))
















