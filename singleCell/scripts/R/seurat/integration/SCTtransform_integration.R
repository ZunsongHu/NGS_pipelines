
options(future.globals.maxSize = 150*1024 * 1024^2) #150GB

# options(future.globals.maxSize = 10*1024 * 1024^2) #10GB

gc(rm(list=ls()))
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")
# setwd("/scratch/zuhu/tmp/GC/RNA_SC")
#temp parameters
# integration_label="BeforeT"
# 
# rds_in="out_sc/EP-01-B-10xT_exp/Seurat/EP-01-B-10xT_exp.SeuratObj.raw.rds,
# out_sc/EP-02-B-10xT_exp/Seurat/EP-02-B-10xT_exp.SeuratObj.raw.rds"
# 
# dir_out=    "out_sc/integration/"

#get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)

integration_label=args[1]
rds_in=args[2]
dir_out=args[3]

n_max=20000

#create output folder -------------------------
if(!dir.exists(dirname(dir_out))){dir.create(dirname(dir_out))}

if(!dir.exists(dir_out)){dir.create(dir_out)}

dir_out_=dir_out

if(!dir.exists(dir_out_)){dir.create(dir_out_)}

rds_out=paste0(dir_out_,"/",integration_label,"_SCT_integrated.rds")

#read obj --------------
message("Reading Obj")
rds_in_list=unlist(strsplit(rds_in,"[,\n]"))
rds_in_list=rds_in_list[!rds_in_list==""]

n_obj=length(rds_in_list)

obj_list=list()

for(i in 1:n_obj){
  message(paste0("Reading ",rds_in_list[i]))
  obj_=readRDS(rds_in_list[i])
  obj_list[[i]]=obj_
}
message("Done Reading Obj\n")


message("subseting\n")
for(i in 1:n_obj){
  set.seed(10)
  barcode_sample=get_barcode_df(obj_list[[i]]) 
  
  if(nrow(barcode_sample)>n_max){
    barcode_sample=barcode_sample %>% sample_n(n_max)
  } else {
    barcode_sample=barcode_sample
  }
  
  obj_list[[i]]=subset(obj_list[[i]],cells=barcode_sample$barcode)
}
message("Done subseting\n")


#SCTransform ------------------------
# print("Performing SCTransform")
# for(i in 1:n_obj){
#   obj_list[[i]] <- SCTransform(obj_list[[i]], verbose = FALSE,variable.features.n = 3000)
# }

#integration --------------------------
# library(future)
# plan("multiprocess", workers = 4)
message("Integrating\n")
features_integration=SelectIntegrationFeatures(object.list = obj_list, nfeatures = 10000)

obj_list = PrepSCTIntegration(object.list = obj_list, anchor.features = features_integration,verbose = FALSE)

anchors_integration = FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
                                             anchor.features = features_integration, verbose = FALSE)

#,reference = 1

obj_integrated = IntegrateData(anchorset = anchors_integration, normalization.method = "SCT",verbose = FALSE)
message("Done integrating\n")

#run PCA UMAP findCluster ----------------
message("running PCA UMAP findCluster")

obj_integrated = RunPCA(obj_integrated, verbose = FALSE)
obj_integrated = RunUMAP(obj_integrated, dims = 1:50,n.neighbors = 30)

obj_integrated = FindNeighbors(obj_integrated, dims = 1:50, verbose = FALSE)

obj_integrated = FindClusters(obj_integrated, verbose = FALSE,resolution = 2)

#draw plot -------------------
DimPlot(obj_integrated,group.by = "group")
ggsave(paste0(dir_out_,"DimPlot_group",".png"),width=10,height = 10,dpi=300)
ggsave(paste0(dir_out_,"DimPlot_group",".pdf"),width=10,height = 10,dpi=300)

DimPlot(obj_integrated,group.by = "seurat_clusters",label = T,repel = T)
ggsave(paste0(dir_out_,"DimPlot_seurat_clusters",".png"),width=10,height = 10,dpi=300)
ggsave(paste0(dir_out_,"DimPlot_seurat_clusters",".pdf"),width=10,height = 10,dpi=300)

#save output ---------------------
message("Saving output\n")
saveRDS(obj_integrated,file=rds_out)































