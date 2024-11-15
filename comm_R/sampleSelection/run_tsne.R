
#rm(list=ls())

source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

library (SummarizedExperiment, quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# obj_rds="2955/obj_2955.rds"
# dataset_label="test"
# dir_out="test/"
# perplexity=10
# variable_n=200
# feature_panel_name="boruta_general_noCorr"


#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=      args[1]
dataset_label=args[2]
perplexity=   as.numeric(args[3])
variable_n=   as.numeric(args[4])
feature_panel_name=args[5]
dir_out=      args[6]



#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
cat("Reading Obj...\n")
obj_=readRDS(obj_rds)

#run tsne ----
obj_tsne=run_tsne(obj_in = obj_,perplexity_one = perplexity,variable_n = variable_n,feature_panel = feature_panel_name)
df_tsne=obj_tsne[["tsne"]]
df_tsne$COH_sample=row.names(df_tsne)

write_tsv(df_tsne,paste0(dir_out,dataset_label,".perplexityN",perplexity,".FeatureN",variable_n,".tsv"))

draw_DimPlot(obj_tsne,"diag",cols=subtypeCol,axis_by=10,reduction="tsne")
ggsave(paste0(dir_out,dataset_label,".perplexityN",perplexity,".FeatureN",variable_n,".png"),width = 12,height = 11)
ggsave(paste0(dir_out,dataset_label,".perplexityN",perplexity,".FeatureN",variable_n,".pdf"),width = 12,height = 11)































