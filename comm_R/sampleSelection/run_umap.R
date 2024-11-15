
#rm(list=ls())

source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

library (SummarizedExperiment, quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")
# 
# obj_rds="2707/sampleSelection/5/get_topMAD/obj.2707.topMAD5000.rds"
# dir_out="test/"
# dataset_label="2707"
# feature_panel="TopMAD5000"
# variable_n=400
# n_neighbors=10

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=            args[1]
dir_out=            args[2]
dataset_label=      args[3]
feature_panel=      args[4]
variable_n=         as.numeric(args[5])
n_neighbors=        as.numeric(args[6])

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
cat("Reading obj...\n")
obj_=readRDS(obj_rds)

df_info_=as.data.frame(colData(obj_$SE))
names(obj_)
#run tsne ----
obj_umap=run_umap(obj_in = obj_,n_neighbors = n_neighbors,variable_n = variable_n,out_label = "umap",feature_panel = feature_panel)

df_umap=obj_umap[["umap"]]
df_umap$COH_sample=row.names(df_umap)


file_out_base=paste0(dir_out,dataset_label,".",feature_panel,".FeatureN",variable_n,".NeighborsN",n_neighbors,".umap")

cat(paste0(file_out_base,".tsv\n"))
write_tsv(df_umap,paste0(file_out_base,".tsv"))

#draw umap -----------------------------------------------------------------------
p1=draw_DimPlot(obj_in = obj_umap,group.by = "diag",cols=subtypeCol,axis_by=5,reduction="umap")+
  theme(panel.grid.major = element_line(color = "grey",linetype = 1),
        panel.grid.minor = element_line(color = "grey",linetype = 1))
p1
ggsave(paste0(file_out_base,".png"),width = 13,height = 11)
ggsave(paste0(file_out_base,".pdf"),width = 13,height = 11)

p2=draw_DimPlot(obj_umap,"diag_raw",cols=subtypeCol,axis_by=5,reduction="umap")+
  theme(panel.grid.major = element_line(color = "grey",linetype = 1),
        panel.grid.minor = element_line(color = "grey",linetype = 1))
p2
ggsave(paste0(file_out_base,"_withPhlike.png"),width = 13,height = 11)
ggsave(paste0(file_out_base,"_withPhlike.pdf"),width = 13,height = 11)

df_pax5Mutation=table(obj_umap$SE$PAX5) %>% as.data.frame() %>% arrange(desc(Freq)) %>% 
  mutate(Var1=as.character(Var1)) %>% filter(!Var1=="0") %>% slice_head(n=15)

#PAX5 Mutation
col_in=unique(subtypeCol)[1:nrow(df_pax5Mutation)]
names(col_in)=df_pax5Mutation$Var1

p3=draw_DimPlot(obj_umap,"PAX5",cols=col_in,axis_by=5,reduction="umap")+
  theme(panel.grid.major = element_line(color = "grey",linetype = 1),
        panel.grid.minor = element_line(color = "grey",linetype = 1))
p3
ggsave(paste0(file_out_base,"_PAX5Mutation.png"),width = 13,height = 11)
ggsave(paste0(file_out_base,"_PAX5Mutation.pdf"),width = 13,height = 11)



























