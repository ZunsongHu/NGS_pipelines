
gc(rm(list=ls()))

source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# 
vst_rds="1974/obj_1974_raw.rds"
dataset_label="1974"
var_batch="library_1"
dir_out="test/"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

vst_rds=      args[1]
dataset_label=args[2]
var_batch=    args[3]
dir_out=      args[4]

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
cat("Read obj...\n")
obj_rna=readRDS(vst_rds)

#get batch freq plot ----
df_sample_batch=as.data.frame(colData(obj_rna$SE))

batch_value=unlist(df_sample_batch[var_batch])

p1=as.data.frame(table(batch_value)) %>% 
  mutate(batch_value=as.character(batch_value)) %>%
  ggplot(aes(x=batch_value, y=Freq)) +
  geom_bar(stat="identity") +
  geom_label(aes(label=Freq))+
  theme(axis.text.x = element_text(angle = 45, hjust=1))

ggsave(paste0(dir_out,"barplot_batchFreq.png"),width = 7,height = 6)

#correct batch ----

obj_rna=correct_batch_var(obj_in = obj_rna,batch_var = var_batch,assay_name_in = "vst")

cat("Saving output...\n")
saveRDS(obj_rna,paste0(dir_out,"obj.",dataset_label,".batch_correction.rds"))














































