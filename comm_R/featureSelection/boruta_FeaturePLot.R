# rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

library(SummarizedExperiment, quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# obj_rds=            "../tsne/2042_middleGroup/obj_2042MiddleGroup.rds"
# label=              "2042_middleGroup"
# boruta_topFeatures= "2042_middleGroup/Boruta_summarise/Boruta_2042_middleGroup_TopList.tsv"
# dir_out=            "test/"


#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=            args[1]
label=              args[2]
boruta_topFeatures= args[3]
dir_out=            args[4]

update_topFeatureTable=T

RankMin_maxValue=2000
# RankMin_maxValue=1000

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read data ----
obj_=readRDS(obj_rds)
df_borutaTopFeatures=read_tsv(boruta_topFeatures)

obj_$SE=obj_$SE[df_borutaTopFeatures$feature_id,]
saveRDS(obj_,paste0(dir_out,"obj_",label,"_topFeatures.rds"))

#get matrix for included features ----
features_in=df_borutaTopFeatures$feature_id[df_borutaTopFeatures$rank_min <= RankMin_maxValue]

df_vst=assays(obj_$SE)[["vst"]]
df_vst=df_vst[features_in,]

df_info=as.data.frame(colData(obj_$SE))

obj_tmp=list(SE=SummarizedExperiment(assays=list(vst=df_vst), colData = df_info))


#draw boxplot ----
feature=features_in[1]

if(!dir.exists(paste0(dir_out,"boxplot/"))){dir.create(paste0(dir_out,"boxplot/"))}

for(feature in features_in){
  draw_BoxPlot(obj_in = obj_tmp,group.by = "diag",features = feature,useGeneName = F)
  ggsave(paste0(dir_out,"boxplot/",gsub("[,+:]","_",feature),".png"),width=8,height=5)
  ggsave(paste0(dir_out,"boxplot/",gsub("[,+:]","_",feature),".pdf"),width=8,height=5)
}















# 
# print("Remove low expression genes")
# obj_tmp=remove_lowexpression_genes(obj_tmp,assay_name_in="vst",cutoff=5)
# 
# print("Get Non correlated genes")
# nonCorr_genes=get_nonCorr_genes(obj_in = obj_tmp,highCorrCutoff = 0.75)
# 
# df_borutaTopFeatures=df_borutaTopFeatures %>%
#   mutate(in_NonCorr_list=ifelse(gene_id %in% nonCorr_genes,"Yes","No"))
# 
# write_tsv(df_borutaTopFeatures,"2042/Boruta_summarise/Boruta_2042_TopList_withNonCorr.tsv")
# 
# 
# table(df_borutaTopFeatures$in_NonCorr_list)
# 
# obj_$boruta_genes=df_borutaTopFeatures$gene_id[df_borutaTopFeatures$in_NonCorr_list=="Yes"]
# perplexity_one=10
# variable_n=150
# dir.create("2042/Boruta_summarise/tsne")
# for (perplexity_one in c(10,20,30)){
#   for (variable_n in c(150,200,300,400,500,600,700,800,900)){
#     obj_=run_tsne(obj_in = obj_,perplexity_one = perplexity_one,dims = ,variable_n = variable_n,out_label = "tsne",feature_panel = "boruta_genes")
#     draw_DimPlot(obj_in = obj_,group.by = "diag",reduction = "tsne")
#     ggsave(paste0("2042/Boruta_summarise/tsne/tsne_",label,"_Boruta_PerplexityN",perplexity_one,"_FeatureN",variable_n,".png"),width = 12,height = 10)
#     ggsave(paste0("2042/Boruta_summarise/tsne/tsne_",label,"_Boruta_PerplexityN",perplexity_one,"_FeatureN",variable_n,".pdf"),width = 12,height = 10)
#   }
# }
# 
# 
# #get Top feature plot ----
# features_top=df_borutaTopFeatures$gene_id[df_borutaTopFeatures$rank_min <=100]
# 
# dir.create("2042/Boruta_summarise/featurePlot")
# 
# obj_=run_tsne(obj_in = obj_,perplexity_one = 30,dims = ,variable_n = 800,out_label = "tsne",feature_panel = "boruta_genes")
# 
# 
# for(feature in features_top){
#   draw_featurePlot(obj_in = obj_,features = feature,reduction = "tsne")
#   
#   gene_label=df_borutaTopFeatures$gene_name[df_borutaTopFeatures$gene_id==feature]
#   minRank=df_borutaTopFeatures$rank_min[df_borutaTopFeatures$gene_id==feature]
#   ggsave(paste0("2042/Boruta_summarise/featurePlot/tsne_",label,"_MinRank",minRank,"_",gene_label,".png"),width = 12,height = 10)
#   ggsave(paste0("2042/Boruta_summarise/featurePlot/tsne_",label,"_MinRank",minRank,"_",gene_label,".pdf"),width = 12,height = 10)
#   
# }
# 
# na=c("MEF2C", "HDAC9", "MEGF10", "CDX2", "NUTM1", "PBX1", "MEIS1", "FLT4", "HLF", "CRLF2", "PAX5")
# 




















