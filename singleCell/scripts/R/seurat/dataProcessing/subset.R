# rm(list=ls())
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

# #temp parameters ----
# setwd("/scratch/zuhu/tmp/SC_CMV/")

rds_in= "out_raw/integration/integration/integration_SCT_integrated.rds"
label=  "integration"
use_featureGene="T"
file_featureGene="id/HSC2Monocyte_gene.list"
subset_label="HSC2Monocyte"
dir_out="test/"
  

#get parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=args[1]
label=args[2]
use_featureGene=args[3]
file_featureGene=args[4]
subset_label=args[5]
dir_out=args[6]

#read rds
cat("Reading obj...\n")
obj_raw=readRDS(rds_in)

if(!dir.exists(dir_out)){dir.create(dir_out)}

#use feature gene ----
if(use_featureGene=="T"){
  featureGenes=unique(readLines(file_featureGene))
  barcode_list1=get_ExpressionBarcode(obj_in = obj_raw,key_genes = featureGenes,KeepTopN = 500)
  
  obj_subset=subset(obj_raw,cells=barcode_list1$barcode)
  df_matrix_count=as.data.frame(as.matrix(GetAssayData(obj_subset,assay = "RNA",slot = "count")))
  
  df_group=get_metadataVar(inObj = obj_subset,var = "group")
  obj_sub=CreateSeuratObject(df_matrix_count)
  
  obj_sub[["group"]]=df_group$group
  
  #run batch correction ----
  cat("Run batch correction...\n")
  obj_sub=batch_correction_seurat(obj_sub,var_batch="group",n_varaible=6000,method="SCT")
  
  obj_sub=RunUMAP(obj_sub,dims = 1:50)
  
  cat("Saving output...\n")
  saveRDS(obj_sub,paste0(dir_out,"obj_",label,"_",subset_label,".rds"))
  
}














