options(future.globals.maxSize = 10*1024 * 1024^2) #10GB

# gc(rm(list=ls()))
options(warn=-1)
suppressPackageStartupMessages({
  library(ggplot2)
  library(Seurat)
  library(SummarizedExperiment)
})
options(warn=0)

source("/net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/singleCell/scripts/R/singleR/run_singleR.R")

#temp parameters
# sc_id="B265_D1"
# 
# rds_in=     paste0("/scratch/zuhu/project/ZhaohuiGu/scBALL_NBC_2022/out_sc/",sc_id,"_Expr5/Seurat/",sc_id,"_Expr5.SeuratObj.QC.rds")
# dir_out=    paste0("/scratch/zuhu/project/ZhaohuiGu/scBALL_NBC_2022/out_sc/",sc_id,"_Expr5/SingleR_annotation_1ref/")
# 
# ref_rds="/home/zgu_labs/bin/R/SingleCell/SingleR/ref/ref_singleR_20celltype.rds"
# sample="B265_D1"
# label="1M_20celltypes"
# anno_label="celltype"
# col_file="/home/zgu_labs/bin/R/SingleCell/SingleR/ref/col_20celltypes.rds"

#get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=args[1]
dir_out=args[2]
sample=args[3]
label=args[4]
ref_rds=args[5]
anno_label=args[6]
col_file=args[7]


#create output folder -------------------------
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read ref ------------------------
ref_=readRDS(ref_rds)
anno_levels=unique(unlist(colData(ref_)[anno_label]))
col_in=readRDS(col_file)

#read seurat obj --------------
obj_QC=readRDS(rds_in)
print(obj_QC)

# set.seed(10)
# df_barcode=get_barcode_df(obj_QC) %>% sample_n(500)
# obj_QC=subset(obj_QC,cells=df_barcode$barcode)

obj_QC=FindNeighbors(obj_QC)
obj_QC=FindClusters(obj_QC,resolution = 2)
DimPlot(obj_QC,group.by = "seurat_clusters")
ggsave(paste0(dir_out,"DimPlot_SeuratCluster",".png"),width=10,height = 10,)

#draw feature plot --------------------------------------
key_features=c("GNLY", "CD3E", "CD8A", "IL7R", "PAX5","CD19","CD34","MS4A1","PBX1","MEIS1","TCF3")

FeaturePlot(obj_QC,features = key_features,ncol = 4)
ggsave(paste0(dir_out,"FeaturePlot",".png"),width=24,height = 12,)

#annotation -------------------------

singleR_out_=run_singleR(obj_QC,ref_,anno_label)

obj_QC=draw_feature_plot_fromDF(obj_QC,singleR_out_$df_singleR_out,'labels',paste0(sample,".",label,".Celllevel"),col_in,dir_out)
obj_QC=draw_feature_plot_fromDF(obj_QC,singleR_out_$df_singleR_out,'labels_cluster',paste0(sample,".",label,".Clusterlevel"),col_in,dir_out)

file_singleR_out=paste0(dir_out,sample,".",label,".SingleR_out.rds")

print(file_singleR_out)

saveRDS(singleR_out_,paste0(file_singleR_out))

#get stringent labels--------------------
# df_label_stringent=get_barcode_df(obj_QC) %>% 
#   left_join(get_embedding_df(obj_QC,"umap")) %>% 
#   left_join(get_stringent_label(singleR_out_$out_singleR,0.75))
# 
# write_tsv(df_label_stringent,paste0(dir_out,"df_label_stringent",".tsv"))
# 
# ggplot(df_label_stringent[!is.na(df_label_stringent$label_stringent),],aes(x=UMAP_1,y=UMAP_2)) + 
#   geom_point(aes(color=label_stringent)) + theme_bw() + ggtitle("_Stringent") +
#   scale_color_manual(values = col_in[names(col_in) %in% df_label_stringent$label_stringent]) 
# 
# ggsave(paste0(dir_out,"SingleR_annotation_","_Stringent",".png"),width=10,height = 10)
# ggsave(paste0(dir_out,"SingleR_annotation_","_Stringent",".pdf"),width=10,height = 10)

#plot delta ----------------
SingleR::plotDeltaDistribution(singleR_out_$out_singleR)
ggsave(paste0(dir_out,sample,".",label,".DeltaDistribution_cellLevel",".png"),width=15,height = 12)

SingleR::plotDeltaDistribution(singleR_out_$out_singleR_cluster)
ggsave(paste0(dir_out,sample,".",label,".DeltaDistribution_ClusterLevel",".png"),width=15,height = 12)

#output  ---------------

saveRDS(obj_QC,file=paste0(dir_out,sample,".",label,".Seurat.rds"))
























