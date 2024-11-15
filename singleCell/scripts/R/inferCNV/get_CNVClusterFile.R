#set ----------------------
gc(rm(list=ls()))
library(infercnv)
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

#COH002922 ----
setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq/")
rds_in="out_sc/COH002922_D1/Seurat/COH002922_D1.SeuratObj.QC.rds"
dir_out="out_sc/COH002922_D1/CNV_manual/"
if(!dir.exists(dir_out)){dir.create(dir_out)}


obj_=readRDS(rds_in)
var_for_analysis="seurat_clusters"

DimPlot(obj_,group.by = var_for_analysis,label=T,repel = T)
ggsave(paste0(dir_out,'DimPlot_',var_for_analysis,".png"),width=7,height = 5)

df_for_CNV=get_barcode_df(obj_) %>% left_join(get_embedding_df(obj_,"umap")) %>%
  left_join(get_metadataVar(obj_,"seurat_clusters")) %>% 
  mutate(
    cluster_for_CNV=ifelse(seurat_clusters %in% c(0,1,2,3),1,ifelse(seurat_clusters==4,2,ifelse(seurat_clusters==5,0,NA)))
  )
  
table(is.na(df_for_CNV$cluster_for_CNV))

obj_[["cluster_for_CNV"]]=df_for_CNV$cluster_for_CNV


var_for_analysis="cluster_for_CNV"

DimPlot(obj_,group.by = var_for_analysis,label=T,repel = T)
ggsave(paste0(dir_out,'DimPlot_',var_for_analysis,".png"),width=7,height = 5)

write_tsv(df_for_CNV,paste0(dir_out,"df_for_CNV.tsv"))




#COH002917 ----
setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq/")
rds_in="out_sc/COH002917_R1/Seurat/COH002917_R1.SeuratObj.QC.rds"
dir_out="out_sc/COH002917_R1/CNV_manual/"
df_CRLF2=read_tsv("out_raw/COHSC000004_Genotyping3/count/count_CB.tsv")
if(!dir.exists(dir_out)){dir.create(dir_out)}

obj_=readRDS(rds_in)
var_for_analysis="seurat_clusters"

DimPlot(obj_,group.by = var_for_analysis,label=T,repel = T)
ggsave(paste0(dir_out,'DimPlot_',var_for_analysis,".png"),width=7,height = 5)

df_for_CNV=get_barcode_df(obj_) %>% left_join(get_embedding_df(obj_,"umap")) %>%
  left_join(get_metadataVar(obj_,"seurat_clusters")) %>% 
  left_join(df_CRLF2 %>% mutate(barcode=V1)) %>%
  mutate(
    seurat_clusters=as.numeric(as.character(seurat_clusters)),
    cluster_for_CNV=ifelse(seurat_clusters %in% c(1,2,6),1,
                           ifelse(seurat_clusters==7,6,seurat_clusters)),
    cluster_for_CNV=ifelse(!cluster_for_CNV==3,cluster_for_CNV,
                           ifelse(count >0,3,
                                  ifelse(count ==0,2,cluster_for_CNV))),
    cluster_for_CNV=ifelse(cluster_for_CNV==2 & UMAP_2 < -1.5,0,cluster_for_CNV),
    cluster_for_CNV=ifelse(cluster_for_CNV==0 & UMAP_2 > -1.5,2,cluster_for_CNV),
    
    cluster_for_CNV=ifelse(cluster_for_CNV==2 & UMAP_2 < 0,2.1,cluster_for_CNV),
    cluster_for_CNV=ifelse(cluster_for_CNV==2 & UMAP_2 < 1,2.2,cluster_for_CNV),
    
    
  ) %>% select(barcode,cluster_for_CNV)


table(is.na(df_for_CNV$cluster_for_CNV))

obj_[["cluster_for_CNV"]]=df_for_CNV$cluster_for_CNV

var_for_analysis="cluster_for_CNV"

DimPlot(obj_,group.by = var_for_analysis,label=T,repel = T)
ggsave(paste0(dir_out,'DimPlot_',var_for_analysis,".png"),width=7,height = 5)

write_tsv(df_for_CNV,paste0(dir_out,"df_for_CNV.tsv"))



#temp parameters
setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq/")
rds_in="out_sc/M0017/Seurat/M0017.SeuratObj.QC.rds"
dir_out="out_sc/M0017/CNV/"
if(!dir.exists(dir_out)){dir.create(dir_out)}

obj=readRDS(rds_in)

table(obj$seurat_clusters)

df_cluster=get_metadataVar(obj,"seurat_clusters")
names(df_cluster)[2]="cluster_for_CNV"

write_tsv(df_cluster,paste0(dir_out,"df_for_CNV.tsv"))


rds_in="out_sc/COH002917_R1/Seurat/COH002917_R1.SeuratObj.QC.rds"
dir_out="out_sc/COH002917_R1/CNV_test/"
var_for_analysis="cluster_CRLF2"
ref_level="0"
species="human"
file_var="out_sc/COH002917_R1/CNV_test/cluster_CRLF2.tsv"


obj_=readRDS(rds_in)

obj_=FindClusters(obj_,resolution = 1)

FeaturePlot(obj_,features = "CRLF2")
ggsave(paste0(dir_out,"umap_CRLF2",".png"),height=10,width = 10)


df_CRLF2=read_tsv("out_raw/COHSC000004_Genotyping3/count/count_CB.tsv")

names(df_CRLF2)[1]="barcode"
ggsave(paste0(dir_out,"umap_CNVClusters",".png"),height=10,width = 10)

df_umap_CRLF2=get_embedding_df(obj_,"umap") %>%
  left_join(df_CRLF2) %>%
  mutate(
    cluster_for_CNV=ifelse(UMAP_1< -10,1,
                           ifelse(UMAP_1< -6,2,
                                  ifelse(UMAP_1< -3.7,3,
                                         ifelse(UMAP_2< -1.5,4,5)))),
    cluster_for_CNV=ifelse(is.na(count),cluster_for_CNV,
                           ifelse(count > 0,6,cluster_for_CNV))
  )

table(df_umap_CRLF2$cluster_for_CNV)/sum(table(df_umap_CRLF2$cluster_for_CNV))

sum(table(df_umap_CRLF2$cluster_for_CNV))

obj_[["cluster_for_CNV"]]=df_umap_CRLF2$cluster_for_CNV

DimPlot(obj_,group.by = "cluster_for_CNV",label=T,repel = T)
ggsave(paste0(dir_out,"umap_CNVClusters",".png"),height=10,width = 10)

write_tsv(df_umap_CRLF2,paste0(dir_out,"df_for_CNV.tsv"))


x=read.table("out_sc/COH002917_R1/CNV/annotations_file.tsv",sep="\t",stringsAsFactors = F)

table(x$V2)

sum(table(df_umap_CRLF2$cluster_for_CNV))

df_umap_CRLF2[is.na(df_umap_CRLF2$cluster_for_CNV),]





