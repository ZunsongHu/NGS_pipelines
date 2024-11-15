
gc(rm(list=ls()))

source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

library (SummarizedExperiment, quietly = T)

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

vst_rds=      args[1]
dataset_label=args[2]
dir_out=      args[3]

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/2863/")
# vst_rds="batch_correction/batch_correction.mini.rds"
# dataset_label="test"
# dir_out="test/"

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read data ----
vst_matrx=readRDS(vst_rds)

feature_all=row.names(vst_matrx)

df_id=data.frame(COH_sample=colnames(vst_matrx),stringsAsFactors = F)
row.names(df_id)=colnames(vst_matrx)

obj_=list(SE=SummarizedExperiment(list(vst=vst_matrx), colData = df_id))

#run mad ----
print("codingNoXYM")
obj_codingNoXYM=obj_
obj_codingNoXYM$SE=obj_codingNoXYM$SE[gene_in,]

print("remove low exp")
obj_codingNoXYM_highE=remove_lowexpression_genes(obj_in = obj_codingNoXYM,cutoff = 5)

print("top MAD10000")
TopMAD10000=get_variable_genes(obj_in=obj_codingNoXYM_highE,N_genes=10000,method="MAD")

obj_codingNoXYM_highE_MAD10000=obj_codingNoXYM_highE
obj_codingNoXYM_highE_MAD10000$SE=obj_codingNoXYM_highE_MAD10000$SE[TopMAD10000,]

print("remove corr")
genes_NoCorr=get_nonCorr_genes(obj_in = obj_codingNoXYM_highE_MAD10000,highCorrCutoff = 0.75)
obj_codingNoXYM_highE_MAD10000_noCorr=obj_codingNoXYM_highE_MAD10000
obj_codingNoXYM_highE_MAD10000_noCorr$SE=obj_codingNoXYM_highE_MAD10000_noCorr$SE[genes_NoCorr,]

print("top MAD5000")
TopMAD5000=get_variable_genes(obj_in=obj_codingNoXYM_highE_MAD10000_noCorr,N_genes=5000,method="MAD")
obj_MAD5000=obj_codingNoXYM_highE_MAD10000_noCorr
obj_MAD5000$SE=obj_MAD5000$SE[TopMAD5000,]

obj_MAD5000[["variable_genes"]]=TopMAD5000

cat(" ")
print(paste0("FeatureNum, SampleNum: ",paste0(dim(assays(obj_MAD5000$SE)$vst),collapse = ", ")))

print("saving RDS")
saveRDS(obj_MAD5000,paste0(dir_out,"obj_MAD5000.",dataset_label,".rds"))



