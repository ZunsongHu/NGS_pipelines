options(future.globals.maxSize = 10*1024 * 1024^2) #2GB

rm(list=ls())
library(ggplot2)
library(Seurat)
library(SummarizedExperiment)
library(cowplot)
library(ggrepel)

source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

# parameters temp ----
# setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq_PublicData/scBALL_NCB_2022/")
# rds_in=       'dir_obj/B328_D1_Expr5.SeuratObj.QC.rds'
# dir_out=      'out_sc/B328_D1_Expr5/DE/'
# id=           "B328_D1_Expr5"
# file_feature= "out_sc/B328_D1_Expr5/DE/B328_D1_Expr5.CRLF2grorup.tsv"

varGroup4DE="CRLF2g"

# parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=     args[1]
dir_out=    args[2]
id=         args[3]
file_feature=args[4]

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#get count matrix ----
cat("Reading single cell object...\n")

obj_=readRDS(rds_in)

# obj_=subset(obj_,cells=row.names(obj_@meta.data %>% sample_n(300)))

#read feature ----
df_feature_in=get_barcode_df(obj_) %>% left_join(read_tsv(file_feature))
df_feature_in=df_feature_in[c("barcode",varGroup4DE)]
names(df_feature_in)[2]='varGroup4DE'

df_feature_in1=df_feature_in %>% filter(!is.na(varGroup4DE))

#get obj subset ----
obj_subset=subset(obj_,cells=df_feature_in1$barcode)
obj_subset$varGroup4DE=(get_barcode_df(obj_subset) %>% left_join(df_feature_in1))$varGroup4DE

#get feature dimplot ----
p1=DimPlot(obj_subset,group.by = "varGroup4DE")
p1
ggsave(paste0(dir_out,"DimPlot.",id,".",varGroup4DE,".png"),width=9,height = 8)
ggsave(paste0(dir_out,"DimPlot.",id,".",varGroup4DE,".pdf"),width=9,height = 8)

#run DE ----
Idents(obj_subset)=obj_subset$varGroup4DE

df_DE=FindAllMarkers(obj_subset, min.pct = 0.25, logfc.threshold = 0.25)
write_tsv(df_DE,paste0(dir_out,"df_DE.",id,".",varGroup4DE,".tsv"))


levels=as.character(sort(unique(df_DE$cluster)))


#get volcano plot ----
for(level in levels){
  df_DE_=df_DE %>% filter(cluster==level)
  df_DE_label=df_DE_ %>% filter(p_val_adj < 0.05)
  
  
  p1=ggplot(data=df_DE_,aes(x=avg_log2FC,y=-log(p_val_adj))) +
    geom_point() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color="grey77"),
          panel.grid.minor.y = element_line(color="grey77"),
          panel.background = element_rect(color="white",fill="white"),
          axis.line = element_line(colour = "black"),)+
    geom_hline(aes(yintercept=-log10(0.05)),linetype=4,color="red")+
    geom_label_repel(data=df_DE_label,aes(x=avg_log2FC,y=-log(p_val_adj),label=gene))+
    geom_point(data=df_DE_label,color="red") +
    labs(x="Fold change",y="-log10(P)",title = paste0("Volcano plot for ",varGroup4DE,": ",level),min.segment.length=1)
  
  p1
  ggsave(paste0(dir_out,"VolcanoPlot.",id,".",varGroup4DE,".",level,".png"),width = 8,height = 7)
  ggsave(paste0(dir_out,"VolcanoPlot.",id,".",varGroup4DE,".",level,".pdf"),width = 8,height = 7)
}



















# ids=c("B328_D1_Expr5","B590_D1_Expr5","B590_R1_Expr5","B734_D1_Expr5","B887_D1_Expr5")
# ids="B590_X1_Expr5"
# 
# setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq/")
# id="COH002917_R1"
# for(id in ids){
  # rds_in="../../SingleCellSeq/out_sc/COH002917_R1/Seurat/COH002917_R1.SeuratObj.QC.rds"
  # obj_=readRDS(rds_in)
  # obj_
  # dir_out=paste0("out_sc/",id,"/DE/")
  # 
  # if(!dir.exists(dir_out)){dir.create(dir_out)}
  # 
  # file_out=paste0("out_sc/",id,"/DE/",id,".CRLF2grorup.tsv")
  # 
  # df_feature=as.data.frame(t(GetAssayData(obj_)))
  # 
  # df_feature$barcode=row.names(df_feature)
  # 
  # df_feature1=df_feature %>% select(barcode,CRLF2) %>% left_join(get_embedding_df(inobj = obj_,reduction_method = "umap"))
  # 
  # df_feature1=df_feature1 %>% filter(UMAP_1 > -4.5 & UMAP_2< -2)
  # 
  # df_feature_high=df_feature1 %>% filter(CRLF2>0) %>% mutate(CRLF2g="1HighCRLF2")
  # set.seed(10)
  # df_feature_low=df_feature1 %>% filter(CRLF2==0) %>% sample_n(nrow(df_feature_high)) %>% mutate(CRLF2g="0LowCRLF2")
  # df_feature_CRLF2=bind_rows(df_feature_high,df_feature_low)
  # 
  # write_tsv(df_feature_CRLF2,file_out)
  
# }















