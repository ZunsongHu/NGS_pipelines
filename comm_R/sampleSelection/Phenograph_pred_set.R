gc(rm(list=ls()))
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

# library(SummarizedExperiment, quietly = T)
# library(caret,quietly = T)

#testing parameters ----
setwd("/scratch/zuhu/project/BALL/tsne/")
file_list=    "2655/sampleSelection/40/summsize/TopMAD5000/PhenographPred_set.2655.TopMAD5000.list"
dataset_label="2655"
dir_out=      "test/"
cutoff_diag_n= 5

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_list=    args[1]
dataset_label=args[2]
dir_out=      args[3]
cutoff_diag_n= as.numeric(args[4])

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read files list ----
files=readLines(file_list)

#get kept level ----
df_one=read_tsv(files[1]) %>% 
  group_by(diag) %>% 
  mutate(freq_diag=n()) %>% 
  select(diag,freq_diag) %>% 
  distinct()

diag_level_keep=df_one$diag[df_one$freq_diag <= cutoff_diag_n]

#get output ----
df_out=bind_rows(lapply(files, function(file_one){
  read_tsv(file_one)
})) %>% arrange(COH_sample)

df_filter_list=df_out %>% filter(!diag == diag_pred) %>% group_by(COH_sample) %>% mutate(
  FeatureN_failed=paste0(FeatureN,collapse = ","),NegihborN_failed=paste0(top_neighborN,collapse = ",")
) %>% select(COH_sample,diag,FeatureN_failed,NegihborN_failed) %>% distinct() %>% arrange(COH_sample) %>% 
  mutate(source="Phenograph")

# df_out1=df_out %>% filter(
#   (!COH_sample %in% df_filter_list$COH_sample) | diag %in% diag_level_keep
# ) %>% select(COH_sample,diag) %>% distinct() %>% arrange(COH_sample)

# table(df_out1$diag)
output_base=paste0(dir_out,"PhenographPred_set.",dataset_label)

write_tsv(df_out,paste0(output_base,".all.tsv"))
write_tsv(df_filter_list,paste0(output_base,".nonMatched.tsv"))
# write_tsv(df_out1,paste0(output_base,".matched.tsv"))
















