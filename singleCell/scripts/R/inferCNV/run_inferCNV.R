#set ----------------------
gc(rm(list=ls()))
library(infercnv)
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

#temp parameters
# setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq/")
# rds_in="out_sc/COH002922_D1/Seurat/COH002922_D1.SeuratObj.QC.rds"
# dir_out="out_sc/COH002922_D1/CNV_test/"
# ref_cluster=11
# species="human"
# var_for_analysis="seurat_clusters"

#get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=args[1]
dir_out=args[2]
var_for_analysis=args[3]
ref_level=args[4]
species=args[5]
file_var=args[6]

#create output folder --------------------------------------------------
if(!dir.exists(dirname(dir_out))){dir.create(dirname(dir_out))}
if(!dir.exists(dir_out)){dir.create(dir_out)}

#file gene list ---------------
if(species=="human"){
  file_geneid="/home/zgu_labs/bin/R/RNAseq/GEP_prediction/Homo_sapiens.GRCh38.V102.summary.txt"
}

if(species=="mouse"){
  file_geneid="/ref_genomes/genome_anno/mouse/gtf/V102/Mus_musculus.GRCm38.V102.withChr.summary.txt"
}

df_geneid=read.table(file_geneid,header = T,stringsAsFactors = F)

#read obj -----
obj_=readRDS(rds_in)

if(!(file_var==""|is.na(file_var))){obj_[[var_for_analysis]]=unlist(read_tsv(file_var)[var_for_analysis])}

# obj_=FindNeighbors(obj_,k.param = 3)
# obj_=FindClusters(obj_,resolution = 2)
# 
# df_umap=get_embedding_df(obj_,"umap") %>%
#   mutate(seurat_clusters=ifelse(UMAP_2> -5,1,ifelse(UMAP_1 < 0,2,3)))
# 
# obj_[["seurat_clusters"]]=df_umap$seurat_clusters

DimPlot(obj_,group.by = var_for_analysis,label=T,repel = T)
ggsave(paste0(dir_out,'DimPlot_',var_for_analysis,".png"),width=7,height = 5)

#get  raw_counts_matrix ------------
counts_=as.matrix(GetAssayData(obj_,slot="counts"))

genes_intersect=intersect(row.names(counts_),df_geneid$gene_name)
df_geneid_in=df_geneid %>% filter(gene_name %in% genes_intersect)
df_gene_freq=as.data.frame(table(df_geneid_in$gene_name))
genes_intersect_in=as.character(df_gene_freq$Var1[df_gene_freq$Freq==1])

counts_in=counts_[genes_intersect_in,]
print(dim(counts_in))

#get gene_order_file ---------------
df_geneid_in=df_geneid[df_geneid$gene_name %in% genes_intersect_in,] 
row.names(df_geneid_in)=df_geneid_in$gene_name
df_geneid_in=df_geneid_in[row.names(counts_in),] %>% 
  mutate(chr1=gsub("chr","",chr),
         chr1=ifelse(chr1 %in% c('1','2','3','4','5','6','7','8','9'),sprintf("%02d",as.numeric(chr1)),chr1)) %>%
  arrange(chr1,start) %>% select(gene_name,chr,start,end) 

write.table(df_geneid_in,paste0(dir_out,'gene_order_file',".tsv"),row.names = F,col.names = F,sep="\t",quote = F)

#get annotations_file ----------
df_cluster=get_metadataVar(obj_,var_for_analysis) 
names(df_cluster)[2]="cluster"
df_cluster$cluster=paste0("cluster",as.character(df_cluster$cluster))

write.table(df_cluster,paste0(dir_out,'annotations_file',".tsv"),col.names = F,row.names = F,sep="\t",quote = F)


infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_,
                                    delim="\t",
                                    annotations_file=paste0(dir_out,'annotations_file',".tsv"),
                                    gene_order_file=paste0(dir_out,'gene_order_file',".tsv"),
                                    ref_group_names=paste0("cluster",ref_level),
                                    chr_exclude=c('chrY', 'chrM')) 

saveRDS(infercnv_obj,file = paste0(dir_out,'infercnv_obj_raw',".rds"))


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=dir_out, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             num_threads=8)

saveRDS(infercnv_obj,file = paste0(dir_out,'infercnv_obj_out',".rds"))

















