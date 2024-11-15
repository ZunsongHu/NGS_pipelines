
rm(list=ls())
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")
setwd("/scratch/zuhu/tmp/SC_CMV/")

dataset_label="YN"
dataset_label="YY"
dataset_label="NY"
obj_rds=paste0("out_raw/",dataset_label,"/Seurat/",dataset_label,".SeuratObj.QC.rds")
dir_out=paste0("out_raw/",dataset_label,"/velocyto_run10X/")

dataset_label="integration"
obj_rds="out_raw/integration/integration/integration_SCT_integrated.rds"
dir_out="out_raw/integration/integration/velocity/"

if(!dir.exists(dir_out)){dir.create(dir_out)}

resolution=0.1

obj_=readRDS(obj_rds)
obj_=FindClusters(obj_,resolution = resolution)

df_umap=get_embedding_df(obj_,"umap") %>% 
  left_join(get_metadataVar(obj_,"seurat_clusters")) %>%
  left_join(get_metadataVar(obj_,"group")) %>%
  mutate(barcode1=ifelse(group=="YN",paste0(barcode,"-0"),
                         ifelse(group=="YY",paste0(barcode,"-1"),
                                ifelse(group=="NY",paste0(barcode,"-2"),NA))))

table(word(df_umap$barcode,2,sep="_"),df_umap$group)

# cols_all=c('gold2','pink','#1F78B5','#FF00B6','#00FF00','magenta3',
#            'red4','#FFA620','#808000','orangered','#5A5156','grey40','#16FF32',
#            'darkgoldenrod4','#A8DD00','#66C2A6','seagreen2','black',
#            'skyblue','#3E9F32','#F6222E','#1E90FF','#E4E1E3','blue3','lightslateblue',
#            '#CCCC33','cyan','grey75','#D27B1C86','#E6BEFF','#469990')

cols_all=c('#0000FF','#FF0000','#00FF00','#000033','#FF00B6','#005300','#FFD300',
           '#009FFF','#9A4D42','#00FFBE','#783FC1','#1F9698','#FFACFD','#B1CC71',
           '#F1085C','#FE8F42','#DD00FF','#201A01','#720055','#766C95','#02AD24',
           '#C8FF00','#886C00','#FFB79F','#858567','#A10300','#14F9FF')

length(unique(df_info$barcode))

cols_in=cols_all[1:length(levels(df_umap$seurat_clusters))]
names(cols_in)=levels(df_umap$seurat_clusters)

df_col=data.frame(seurat_clusters=levels(df_umap$seurat_clusters),colour=cols_in,stringsAsFactors = F)

DimPlot(obj_,group.by = "seurat_clusters",cols = cols_in,label = T,repel = T)
ggsave(paste0(dir_out,dataset_label,"_seurat_clusters.png"),width = 12,height = 10)
ggsave(paste0(dir_out,dataset_label,"_seurat_clusters.pdf"),width = 12,height = 10)


df_info=df_umap %>% left_join(df_col)


write_tsv(df_info,paste0(dir_out,dataset_label,"_info.tsv"))























