# gc(rm(list=ls()))
# rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

library (SummarizedExperiment, quietly = T)
library(UBL, quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# 
# obj_rds=          "../tsne/2042_mRNA/obj_2042mRNA.rds"
# dataset_label =   "2042_mRNA"
# var_group_smote=  "diag"
# smote_N=          50
# dir_out=          "test/"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=          args[1]
dataset_label=    args[2]
var_group_smote=  args[3]
smote_N=          as.numeric(args[4])
dir_out=          args[5]



features=NULL

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
cat("Reading Obj\n")
obj_=readRDS(obj_rds)


#get features ----
cat("Get feature names\n")
if(is.null(features)){features_in=rownames(obj_$SE)}
if(!is.null(features)){features_in=features}
# features_in=features_in[1:100]

#get df for smote ----
cat("Get feature matrix\n")
value_group=unlist(colData(obj_$SE)[var_group_smote])
table(value_group)
obj_$SE[["group"]]=as.factor(value_group)
df_feature_group=get_features_df(obj_,assay_name_in="vst",features=c("group",features_in))

#run smote ----
cat("Run smote\n")
subtypeFreq=table(value_group)
subtypeFreqRatio=as.list(smote_N/subtypeFreq)

set.seed(10)
exprModelSmoted = SmoteClassif(group ~ ., df_feature_group, C.perc = subtypeFreqRatio,k=3)

names(exprModelSmoted)[match("group",names(exprModelSmoted))]=var_group_smote

cat("Save to output\n")
saveRDS(exprModelSmoted,paste0(dir_out,"smote_",dataset_label,"_",var_group_smote,"_EachGroupN",smote_N,".rds"))


#get tsne plot ----

# for (fileone in filelist){
#   print(fileone)
#   df_vst=readRDS(fileone)
#   
#   df_coldata=df_vst["diag"]
#   
#   matrix_=t(df_vst[gene_in])
#   
#   obj_=list(SE=SummarizedExperiment(assays=list(vst=matrix_), colData = df_coldata))
#   
#   print("Remove low expression genes")
#   obj_=remove_lowexpression_genes(obj_,assay_name_in="vst",cutoff=5)
#   
#   print("Get Non correlated genes")
#   nonCorr_genes=get_nonCorr_genes(obj_in = obj_,highCorrCutoff = 0.75)
#   obj_$SE=obj_$SE[nonCorr_genes,]
#   
#   print("Get top MAD genes")
#   variable_genes=get_variable_genes(obj_in = obj_,N_genes = 5000,method = "mad")
#   obj_[["variable_genes"]]=variable_genes
#   
#   obj_=run_tsne(obj_in = obj_,perplexity_one = 30,dims = 2,variable_n = 800,out_label = "tsne")
#   
#   draw_DimPlot(obj_,group.by="diag",cols=subtypeCol,axis_by=10,reduction="tsne")
#   ggsave(paste0("2042/smote/",gsub(".rds","",basename(fileone)),".tsne",".png"),width = 12,height = 10)
#   ggsave(paste0("2042/smote/",gsub(".rds","",basename(fileone)),".tsne",".pdf"),width = 12,height = 10)
# }























