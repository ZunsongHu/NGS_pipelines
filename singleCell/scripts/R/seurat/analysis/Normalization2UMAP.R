# rm(list=ls())
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")


# parameters temp ----
setwd("/scratch/zuhu/tmp/SC_CMV/")

rds_in="out_raw/integration/integration/subset/HSC2Monocyte/obj_integration_HSC2Monocyte.rds"
method="SCT"
method="ordinary"
batch_var="group"
batch_var="NULL"

rds_out="test/test.rds"

#read obj ----
cat("Reading obj...\n")
obj_raw=readRDS(rds_in)

if(!dir.exists(dirname(rds_out))){dir.create(dirname(rds_out))}

if(method=="SCT" & batch_var=="NULL"){
  obj_raw=SCTransform(obj_raw)
}

if(method=="ordinary" & batch_var=="NULL"){
  obj_raw=ScaleData(obj_raw)
  obj_raw=NormalizeData(obj_raw)
  obj_raw=FindVariableFeatures(obj_raw)
}


if(method=="SCT" & !batch_var=="NULL"){
  obj_raw <- SCTransform(
    obj_raw,
    vars.to.regress = c(batch_var),
    verbose=TRUE
  )
  
  obj_list=SplitObject(obj_raw,split.by = batch_var)
  obj_list=lapply(obj_list,SCTransform)
  
  features = SelectIntegrationFeatures(object.list = obj_list, nfeatures = 15000)
  
  obj_list = PrepSCTIntegration(object.list = obj_list, anchor.features = features)
  anchors = FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
                                           anchor.features = features)
  obj_nrom <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
}

ref_singleR=readRDS("/home/zgu_labs/bin/R/SingleCell/SingleR/ref/HPCA_generausage_27celltyps.rds")

df_celltypist=read.csv("out_raw/integration/integration/subset/HSC2Monocyte/celltypist/celltypist_raw/",
                       stringsAsFactors = F)

names(df_celltypist)[1]="barcode"

obj_nrom=RunPCA(obj_nrom)
obj_nrom=RunUMAP(obj_nrom,dims = 1:50)
obj_nrom=FindNeighbors(obj_nrom)
obj_nrom=FindClusters(obj_nrom,resolution = 2)
obj_singleR_out=run_singleR(obj_nrom,ref_singleR = ref_singleR,var_group = "celltype")

obj_nrom[["singleR"]]=obj_singleR_out$df_singleR_out$labels

df_celltypist1=get_barcode_df(obj_nrom) %>% left_join(df_celltypist)

obj_nrom[["celltypist"]]=df_celltypist1$predicted_labels

obj_SCT=obj_

obj_SCTbatch=obj_

obj_SCTbatch1=obj_nrom

DimPlot(obj_SCTbatch,group.by = "group")

run_cellCycle=T















