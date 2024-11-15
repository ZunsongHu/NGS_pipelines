# rm(list=ls())
library (SummarizedExperiment, quietly = T)
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/tsne/")

# file_rds="2707/sampleSelection/1/batch_correction/obj.2707.batch_correction.rds"
# dataset_label="2707"
# keep_XYM="F"
# dir_out="test/"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_rds=       args[1]
dataset_label=  args[2]
keep_XYM=       args[3]
dir_out=        args[4]

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read data ----
cat("Reading Obj...\n")
obj_rna=readRDS(file_rds)

# random
# set.seed(10)
# df_info=as.data.frame(colData(obj_rna$SE)) %>% sample_n(250)
# obj_rna$SE=obj_rna$SE[,df_info$COH_sample]
# names(obj_rna)

#subset ----
if(!keep_XYM=="T"){
  print("codingNoXYM")
  obj_codingNoXYM=obj_rna
  obj_codingNoXYM$SE=obj_codingNoXYM$SE[gene_in,]
}

if(keep_XYM=="T"){
  obj_codingNoXYM=obj_rna
}

cat("Remove low exp...\n")
obj_codingNoXYM_highE=remove_lowexpression_genes(obj_codingNoXYM,cutoff = 5)

cat("Get MAD10000...\n")
obj_codingNoXYM_highE$TopMAD10000=get_variable_genes(obj_in =obj_codingNoXYM_highE,N_genes=10000,method="MAD")

cat("Remove corr...\n")
obj_codingNoXYM_highE$TopMAD10000_noCorr=get_nonCorr_genes(obj_in = obj_codingNoXYM_highE,feature_panel = "TopMAD10000",highCorrCutoff = 0.75)

obj_codingNoXYM_highE_top1000MADnoCorr=obj_codingNoXYM_highE
obj_codingNoXYM_highE_top1000MADnoCorr$SE=obj_codingNoXYM_highE_top1000MADnoCorr$SE[obj_codingNoXYM_highE$TopMAD10000_noCorr,]

cat("Get MAD5000...\n")
obj_codingNoXYM_highE_top1000MADnoCorr$TopMAD5000=get_variable_genes(obj_in=obj_codingNoXYM_highE_top1000MADnoCorr,N_genes=5000,method="MAD")

cat("Get output obj...\n")
obj_topMAD5000=obj_rna
obj_topMAD5000$SE=obj_rna$SE[sort(unique(c(obj_codingNoXYM_highE_top1000MADnoCorr$TopMAD5000,obj_rna$genes_1366))),]
obj_topMAD5000$TopMAD5000=obj_codingNoXYM_highE_top1000MADnoCorr$TopMAD5000

cat(paste0("FeatureNum, SampleNum: ",paste0(dim(assays(obj_topMAD5000$SE)$vst),collapse = ", "),"\n"))

cat("saving output...\n")
saveRDS(obj_topMAD5000,paste0(dir_out,"obj.",dataset_label,".topMAD5000.rds"))

#write df_info
write_tsv(as.data.frame(colData(obj_topMAD5000$SE)),paste0(dir_out,"df_info.",dataset_label,".tsv"))












