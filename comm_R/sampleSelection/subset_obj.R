rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

# library(SummarizedExperiment, quietly = T)
# library(caret,quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# 
# rds_obj="2707/obj_2707_beforeBatch.rds"
# dir_out="test/"
# label="2707"
# file_list="2707/sampleSelection/41/selected/nonMatched.list"
# info_base="2707/sampleSelection/41/get_topMAD/df_info.2707.tsv"
# cutoff_diag_n=100

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_obj=      args[1]
dir_out=      args[2]
label=        args[3]
file_list=    args[4]
info_base=    args[5]
cutoff_diag_n=  as.numeric(args[6])

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
cat("Reading obj...\n")
obj_=readRDS(rds_obj)

df_info_base=read_tsv(info_base)

#get kept level ----

df_one=df_info_base %>% 
  group_by(diag) %>% 
  mutate(freq_diag=n()) %>% 
  select(diag,freq_diag) %>% 
  distinct()

diag_level_keep=df_one$diag[df_one$freq_diag <= cutoff_diag_n]


#read df for matching samples ----
files=readLines(file_list)

files=files[!files=="NULL"]

df_NonMatch_= bind_rows(
  lapply(files, function(x){
    df=read_tsv(x) %>% 
      mutate(
        COH_sample=as.character(COH_sample),
        diag=as.character(diag),
        source=as.character(source),
        FeatureN_failed=as.character(FeatureN_failed),
        NegihborN_failed=as.character(NegihborN_failed),
        dataset_i=match(x,files))
  })) %>% select(COH_sample,diag,source,dataset_i) %>% distinct() %>% 
  group_by(COH_sample) %>% 
  mutate(count_=n()) %>% 
  filter(count_>=0) %>% 
  select(COH_sample,diag) %>% distinct() %>% ungroup() %>% 
  filter(!diag %in% diag_level_keep)

write_tsv(df_NonMatch_,paste0(dir_out,"df_NonMatch_.tsv"))

#draw plot for freqs of diag for matching samples ----
df_match_=df_info_base %>% filter(!COH_sample %in% df_NonMatch_$COH_sample) 
p1=df_match_ %>% 
  group_by(diag) %>% 
  mutate(Freq=n()) %>% select(diag,Freq) %>% 
  ggplot(aes(x=diag, y=Freq)) +
  geom_bar(stat="identity") +
  geom_label(aes(x=diag,y=Freq,label=Freq))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
p1
ggsave(paste0(dir_out,"DiagFreqPlot.",label,".selected.png"),width = 7,height = 6)

#get output obj ----
obj_out=obj_

obj_out$SE=obj_out$SE[,df_match_$COH_sample]

#save output ----
file_out=paste0(dir_out,"obj.",label,".selected.rds")

cat(paste0(file_out,"\n"))
saveRDS(obj_out,file_out)










