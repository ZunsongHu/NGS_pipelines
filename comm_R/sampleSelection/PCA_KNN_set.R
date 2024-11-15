#rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

# library(SummarizedExperiment, quietly = T)
# library(caret,quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")

# file_KNN_list=    "2655/sampleSelection/1/summsize/TopMAD5000/PCAKNN_set.2655.TopMAD5000.list"
# dataset_label="test"
# dir_out=      "test/"
# cutoff_diag_n=5

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

dir_out=          args[1]
dataset_label=    args[2]
file_KNN_list=    args[3]
cutoff_diag_n=    as.numeric(args[4])

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

df_filter_list=df_out %>% filter(!diag == KNN_pred) %>% group_by(COH_sample) %>% mutate(
  FeatureN_failed=paste0(FeatureN,collapse = ","),PC_N_failed=paste0(PC_N,collapse = ",")
) %>% select(COH_sample,diag,FeatureN_failed,PC_N_failed) %>% distinct() %>% arrange(COH_sample) %>% 
  mutate(source="PCA")

# df_out1=df_out %>% filter(
#   (!COH_sample %in% df_filter_list$COH_sample) | diag %in% diag_level_keep
# ) %>% select(COH_sample,diag) %>% distinct()

# table(df_out1$diag)
output_base=paste0(dir_out,"PCAKNN_set.",dataset_label)


write_tsv(df_out,paste0(output_base,".all.tsv"))
write_tsv(df_filter_list,paste0(output_base,".nonMatched.tsv"))
# write_tsv(df_out1,paste0(output_base,".matched.tsv"))


















