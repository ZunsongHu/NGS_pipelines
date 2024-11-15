
options(future.globals.maxSize = 100*1024 * 1024^2) #50GB

gc(rm(list=ls()))
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")
library(stringr)
#temp parameters
rds_in="/scratch/zuhu/project/ZhaohuiGu/singlecell/PublicData/1M/out/integration_sample1/SeuratObj.Integration.DoubleAnno.rds,
/scratch/zuhu/project/ZhaohuiGu/singlecell/PublicData/1M/out/integration_sample2/SeuratObj.Integration.DoubleAnno.rds,
/scratch/zuhu/project/ZhaohuiGu/singlecell/PublicData/1M/out/integration_sample3/SeuratObj.Integration.DoubleAnno.rds,
/scratch/zuhu/project/ZhaohuiGu/singlecell/PublicData/1M/out/integration_sample4/SeuratObj.Integration.DoubleAnno.rds,
/scratch/zuhu/project/ZhaohuiGu/singlecell/PublicData/1M/out/integration_sample5/SeuratObj.Integration.DoubleAnno.rds,
/scratch/zuhu/project/ZhaohuiGu/singlecell/PublicData/1M/out/integration_sample6/SeuratObj.Integration.DoubleAnno.rds,
/scratch/zuhu/project/ZhaohuiGu/singlecell/PublicData/1M/out/integration_sample7/SeuratObj.Integration.DoubleAnno.rds,
/scratch/zuhu/project/ZhaohuiGu/singlecell/PublicData/1M/out/integration_sample8/SeuratObj.Integration.DoubleAnno.rds"

n_in_group=50

#get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=args[1]
n_in_group=as.numeric(args[2])

#make out file name ------------
dir_out=    paste0("/scratch/zuhu/project/ZhaohuiGu/singlecell/PublicData/1M/integrated_all/randomsample",n_in_group,"/")
rds_out=paste0(dir_out,"Integrated_samplenumber",n_in_group,".rds")

#create output folder -------------------------
if(!dir.exists(dirname(dir_out))){dir.create(dirname(dir_out))}
if(!dir.exists(dir_out)){dir.create(dir_out)}


#read obj --------------
print("Reading Obj")
rds_in_list=unlist(strsplit(rds_in,"[,\n]+"))
rds_in_list=rds_in_list[!rds_in_list==""]
print(rds_in_list)

n_obj=length(rds_in_list)

obj_list=list()

for(i in 1:n_obj){
  obj_=readRDS(rds_in_list[i])
  obj_list[[i]]=obj_$obj_QC
}

for(i in 1:n_obj){
  set.seed(10)
  barcode_sample=get_metadataVar(obj_list[[i]],"HPCA_Blineage_cellLevel") %>%  
    mutate(random_number=runif(dim(obj_list[[i]])[2])) %>%
    group_by(HPCA_Blineage_cellLevel) %>%
    arrange(random_number) %>%
    slice_head(n=n_in_group)
  obj_list[[i]]=subset(obj_list[[i]],cells=barcode_sample$barcode)
}

#SCTransform ------------------------
print("Performing SCTransform")
for(i in 1:n_obj){
  obj_list[[i]] <- SCTransform(obj_list[[i]], verbose = FALSE,variable.features.n = 3000)
}

#integration --------------------------
# library(future)
# plan("multiprocess", workers = 4)

features_integration=SelectIntegrationFeatures(object.list = obj_list, nfeatures = 15000)

obj_list = PrepSCTIntegration(object.list = obj_list, anchor.features = features_integration,verbose = FALSE)

anchors_integration = FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",anchor.features = features_integration, verbose = FALSE)

obj_integrated = IntegrateData(anchorset = anchors_integration, normalization.method = "SCT",verbose = FALSE)

#run PCA UMAP findCluster ----------------
print("running PCA UMAP findCluster")

obj_integrated = RunPCA(obj_integrated, verbose = FALSE)
obj_integrated = RunUMAP(obj_integrated, dims = 1:50,n.neighbors = 5,)

obj_integrated = FindNeighbors(obj_integrated, dims = 1:50, verbose = FALSE)

obj_integrated = FindClusters(obj_integrated, verbose = FALSE,resolution = 2)

#draw plot -------------------
DimPlot(obj_integrated,group.by = "group")
ggsave(paste0(dir_out,"DimPlot_group",".png"),width=25,height = 10,dpi=300)
ggsave(paste0(dir_out,"DimPlot_group",".pdf"),width=25,height = 10,dpi=300)

obj_integrated[["sample"]]=word(obj_integrated$group,1,sep="_")
DimPlot(obj_integrated,group.by = "sample")
ggsave(paste0(dir_out,"DimPlot_sample",".png"),width=12,height = 10,dpi=300)
ggsave(paste0(dir_out,"DimPlot_sample",".pdf"),width=12,height = 10,dpi=300)

DimPlot(obj_integrated,group.by = "seurat_clusters",label = T,repel = T)
ggsave(paste0(dir_out,"DimPlot_seurat_clusters",".png"),width=10,height = 10,dpi=300)
ggsave(paste0(dir_out,"DimPlot_seurat_clusters",".pdf"),width=10,height = 10,dpi=300)


DimPlot(obj_integrated,group.by = "HPCA_Blineage_cellLevel",label = T,repel = T,
        cols=DiscretePalette(length(unique(obj_integrated$HPCA_Blineage_cellLevel))))
ggsave(paste0(dir_out,"DimPlot_celltype_celllevel_sample",n_in_group,".png"),width=10,height = 10,dpi=300)
ggsave(paste0(dir_out,"DimPlot_celltype_celllevel_sample",n_in_group,".pdf"),width=10,height = 10,dpi=300)


#save output ---------------------
print("save output")
saveRDS(obj_integrated,file=rds_out)






























