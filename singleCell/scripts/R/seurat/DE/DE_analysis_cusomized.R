gc(rm(list=ls()))
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

#temp parameters ----
setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq/")

file_obj="out_sc/M0021/Seurat/M0021.SeuratObj.QC.rds"
file_VDJ="out_sc/M0021/trust4fq_plot/M0021_trust4_cleaned.tsv"
file_mut="out_sc/M0021/Mutation_plot/M0021_Pax5.P80R_umap.mutation.tsv"
file_celltype="out_sc/M0021/SingleR_annotation_1Ref/M0021.mouse_ImmGenData.SingleR_out.rds"
dir_out="out_sc/M0021/DE_customized/"

x_min=-15
x_max=-2.5

y_min=-1
y_max=5

key_level_vdj="NoVDJ"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_obj=args[1]
file_VDJ=args[2]

file_mut=args[3]
file_celltype=args[4]

dir_out=args[5]

x_min=as.numeric(args[6])
x_max=as.numeric(args[7])
y_min=as.numeric(args[8])
y_max=as.numeric(args[9])

key_level_vdj=args[10]


#read data ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

obj_=readRDS(file_obj)
df_vdj=read_tsv(file_VDJ) %>%
  mutate(
    VDJHC1=ifelse(VDJHC1=="",NA,VDJHC1),
    VDJHC2=ifelse(VDJHC2=="",NA,VDJHC2),
    VDJH_VDJL1=ifelse(VDJH_VDJL1=="",NA,VDJH_VDJL1),
    VDJH_VDJL2=ifelse(VDJH_VDJL2=="",NA,VDJH_VDJL2),
    VDJC=ifelse(VDJC=="",NA,VDJC),
    VDJH_VDJL=ifelse(VDJH_VDJL=="",NA,VDJH_VDJL)
  ) %>%
  mutate(
    VDJC=ifelse(is.na(VDJC),"NoVDJ",VDJC),
    VDJH_VDJL=ifelse(is.na(VDJH_VDJL),"NoVDJ",VDJH_VDJL)
  ) %>%
  select(barcode,VDJHC1,VDJHC2,VDJH_VDJL1,VDJH_VDJL2,VDJC,VDJH_VDJL) 

# df_mut=read_tsv(file_mut) %>%
#   select(barcode,gene_mut,total_CB,alt_count,ref_count,MAF,hasMutation)

df_celltype=readRDS(file_celltype)$df_singleR_out %>% mutate(seurat_clusters=NULL)


#subset obj ----

df_umap=get_embedding_df(obj_,"umap") %>%
  mutate(inx=ifelse(UMAP_1 > x_min & UMAP_1 < x_max & UMAP_2 > y_min &  UMAP_2< y_max ,"in","no")
  ) %>% 
  left_join(df_vdj) %>% 
  filter(inx=="in")

obj_subset=subset(obj_,cells=df_umap$barcode[df_umap$inx=="in"])

obj_subset[["VDJH_VDJL"]]=df_umap$VDJH_VDJL

DimPlot(obj_subset,group.by = "seurat_clusters",label = T,repel = T)
ggsave(paste0(dir_out,"umap_seurate_clusters",".png"),width = 12,height = 5)

#run DE ----
Idents(obj_subset)="VDJH_VDJL"

table(df_umap$VDJH_VDJL)

out_de <- FindMarkers(obj_subset, ident.1 = key_level_vdj, ident.2 = "OtherVDJ")

df_DE=as.data.frame(out_de) %>% mutate(gene=row.names(out_de))

write_tsv(df_DE,paste0(dir_out,"df_DE_VDJ",".tsv"))

#save plot ----
head(as.data.frame(df_DE))

gene_top=(df_DE %>% arrange(p_val_adj) %>% slice_head(n=4))$gene

DimPlot(obj_subset,group.by = "VDJH_VDJL") + FeaturePlot(obj_subset,features=gene_top)
ggsave(paste0(dir_out,"VDJ_DEgenes",".png"),width = 22,height = 15)

saveRDS(obj_subset,paste0(dir_out,"obj_subset",".rds"))






