# options(future.globals.maxSize =10*1024 * 1024^2) #10GB
# gc(rm(list=ls()))
# rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

# suppressWarnings(suppressMessages(require(SummarizedExperiment)))

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_errorRate=            args[1]

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# 
# file_errorRate="2042/GA_selection/1_sum/GA_EachGroupN100_NTrain50_NFeatures800_NKNN5_DimsPca50.tsv"

#create dir out ----
dir_out=paste0(dirname(file_errorRate),"/")
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read df error ----
file_out_prefix=word(basename(file_errorRate),1,sep = "[.]")

strs=unlist(strsplit(file_out_prefix,"_"))
n_for_features=as.numeric(gsub("NFeatures","",strs[grepl("NFeatures",strs)]))

df_error=read_tsv(file_errorRate)

ggplot(df_error,aes(x=error_rate_inTestData)) +
  geom_histogram(binwidth=0.001) +
  geom_density()

ggsave(paste0(dir_out,file_out_prefix,"_ModelError_Histogram.png"),height=8,width=8)
ggsave(paste0(dir_out,file_out_prefix,"_ModelError_Histogram.pdf"),height=8,width=8)

#get gene error rate ----
df_error1=df_error %>%
  select(error_rate_inTestData,features_in) %>%
  tidyr::separate_rows(features_in,sep=",") %>%
  group_by(features_in) %>%
  mutate(error_rate_inTestData=mean(error_rate_inTestData)) %>%
  distinct() %>%
  mutate(gene_id=features_in) %>% left_join(info_gtf) %>%
  arrange(error_rate_inTestData)

write_tsv(df_error1,paste0(dir_out,file_out_prefix,"_GeneError.tsv"))

ggplot(df_error1,aes(x=error_rate_inTestData)) +
  geom_histogram(binwidth=0.0001)+
  geom_density()

ggsave(paste0(dir_out,file_out_prefix,"_GeneError_Histogram.png"),height=8,width=8)
ggsave(paste0(dir_out,file_out_prefix,"_GeneError_Histogram.pdf"),height=8,width=8)


#get top df ----
df_error_top=df_error %>% arrange(error_rate_inTestData) %>% slice_head(n=2) %>%
  mutate(model=paste0("Model",1:n())) %>%
  select(features_in,model) %>%
  tidyr::separate_rows(features_in,sep=",") %>% left_join(df_error1)

df_error_top1=df_error_top %>% 
  mutate(randomNum=runif(n())) %>%
  group_by(model) %>%
  arrange(randomNum) %>%
  slice_head(n=n_for_features/2) %>%
  ungroup() %>%
  select(features_in,gene_name) %>%
  distinct()

df_error_top2=df_error_top  %>% 
  mutate(model=NULL) %>% distinct() %>%
  arrange(error_rate_inTestData,gene_name)

write.table(df_error_top1["features_in"],paste0(dir_out,file_out_prefix,"_top.list"),col.names = F,row.names = F,quote = F)

write_tsv(df_error_top2,paste0(dir_out,file_out_prefix,"_topGene.tsv"))
















































