
options(future.globals.maxSize = 150*1024 * 1024^2) #150GB

# options(future.globals.maxSize = 10*1024 * 1024^2) #10GB

rm(list=ls())
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

# temp parameters
setwd("/scratch/zuhu/tmp/SC_CMV/")
integration_label="threesample"

rds_in="out_raw/NY/Seurat/NY.SeuratObj.QC.rds,out_raw/YN/Seurat/YN.SeuratObj.QC.rds,out_raw/YY/Seurat/YY.SeuratObj.QC.rds"

dir_out=    "out_raw/integration/"

#get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)

integration_label=args[1]
rds_in=args[2]
dir_out=args[3]
method="ordinary"


#create output folder -------------------------
if(!dir.exists(dirname(dir_out))){dir.create(dirname(dir_out))}
if(!dir.exists(dir_out)){dir.create(dir_out)}

rds_out=paste0(dir_out,"/",integration_label,"_SCT_integrated.rds")

#read obj --------------
cat("Reading Obj...\n")
rds_in_list=unlist(strsplit(rds_in,"[,\n]"))
rds_in_list=rds_in_list[!rds_in_list==""]

n_obj=length(rds_in_list)

obj_list=list()

for(i in 1:n_obj){
  obj_=readRDS(rds_in_list[i])
  obj_list[[i]]=obj_
  
}

#normalization ----
if(method=="ordinary"){
  
  cat("Runing preparation...\n")
  
  obj_list1=lapply(obj_list, function(x){
    matrix_count=GetAssayData(x[["RNA"]], slot = "counts")
    obj_=CreateSeuratObject(matrix_count,meta.data = x@meta.data)
    obj_=NormalizeData(obj_)
    obj_=FindVariableFeatures(obj_, selection.method = "vst", nfeatures = 10000,verbose = F)
  })
  
  features <- SelectIntegrationFeatures(object.list = obj_list1,nfeatures = 10000)
  
  obj_list2 = lapply(X = obj_list1, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  
  cat("Finding anchors, using rPCA\n")
  obj_anchors = FindIntegrationAnchors(object.list = obj_list2, anchor.features = features, reduction = "rpca",verbose = F)
  
  obj_combined= IntegrateData(anchorset = obj_anchors)
  
  
  DefaultAssay(obj_combined) = "integrated"
  
  obj_combined=ScaleData(obj_combined)
  
  obj_combined=FindVariableFeatures(obj_combined,nfeatures = 2000)
  
  obj_combined=RunPCA(obj_combined)
  
  obj_combined=FindNeighbors(obj_combined)
  
  obj_combined=FindClusters(obj_combined,resolution = 2)
  
  obj_combined=RunUMAP(obj_combined,dims = 1:50)
  
  # saveRDS(obj_combined,"out_raw/integration/threesample/threesample_integrated.rds")
  
  

  obj_=ScaleData(obj_)
  
  
  in_ids=unlist(c(lapply(obj_list1, function(x){unique(x$group)})))
  
  obj_merge=merge(obj_list1, add.cell.ids = in_ids, project = integration_label)
  
}

obj_list1[[2:3]]



pbmc.big <- merge(pbmc3k, y = c(pbmc4k, pbmc8k), add.cell.ids = in_ids, project = "PBMC15K")
pbmc.big


str(obj_@assays$RNA)

matrix_count=as.matrix(GetAssayData(obj_@, slot = "counts"))
dim(matrix_count)










