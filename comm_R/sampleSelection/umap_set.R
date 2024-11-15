
#rm(list=ls())

source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

library (SummarizedExperiment, quietly = T)

# testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# 
# file_rds_obj="2042/obj_2042.rds"
# file_list="2042/umap/coding/2042.tsv.list"
# dataset_label="test"
# feature_panel_name="coding"
# dir_out="test/"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_rds_obj=       args[1]
file_list=          args[2]
dataset_label=      args[3]
feature_panel_name= args[4]
dir_out=            args[5]

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
cat("Read Obj...\n")
obj_=readRDS(file_rds_obj)

df_info=as.data.frame(colData(obj_$SE))

#read data ----
file_list=readLines(file_list)


df_plot=bind_rows(lapply(file_list, read_tsv)) %>% 
  left_join(df_info) %>% distinct()


gg_dimPlot(df_in = df_plot,x = names(df_plot)[1],y = names(df_plot)[2],var_col = "diag",cols = subtypeCol,plot_title = "Separated by FeatureN",split.by = "FeatureN",
           size = 0.3)

width_base=ceiling(sqrt(length(unique(df_plot$n_neighbors))*length(unique(df_plot$FeatureN))))

width_=width_base*10
height_=ceiling(width_/4*3)

ggsave(paste0(dir_out,"umap_merge.",dataset_label,".",feature_panel_name,".pdf"),width = width_,height = height_)
























