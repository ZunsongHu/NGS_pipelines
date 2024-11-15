# rm(list=ls())

library(stringr)
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

#temp parameters ----
# setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq/scBALL_NCB_2022/")
# 
# file_obj="out_sc/H1_G1_Expr5/Seurat/H1_G1_Expr5.SeuratObj.QC.rds"
# file_region="out_sc/H1_G1_Expr5/count_region/H1_G1_Expr5.AL021978.1.tsv"
# dir_out="test/"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_obj=args[1]
file_region=args[2]
dir_out=args[3]

id=gsub(".SeuratObj.QC.rds","",basename(file_obj))

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
obj_=readRDS(file_obj)

#read count region ----
df_count=read_tsv(file_region)

region=names(df_count)[[2]]

#seurat_clusters and gene matrix ----
df_out=get_barcode_df(obj_) %>%
  left_join(df_count) %>% 
  left_join(get_embedding_df(obj_,"umap")) 

# gg_featurePlot(df_in = df_out,x="UMAP_1",y="UMAP_2",var_col = region,cols = c("purple","tan2","red"),size = 0.3)

table(colnames(obj_)==df_out$barcode)

obj_$feature_count=unlist(df_out[2])

FeaturePlot(object = obj_,features = "feature_count",cols = c("purple","tan2","red"),order = T)+
  labs(title = region)
ggsave(paste0(dir_out,id,"_",region,".pdf"),width=7,height = 6)
ggsave(paste0(dir_out,id,"_",region,".png"),width=7,height = 6)

#save output file -----------------
write_tsv(df_out,paste0(dir_out,"df_regionCount_",id,".",region,".tsv"))





