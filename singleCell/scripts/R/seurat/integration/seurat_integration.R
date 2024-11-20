# usage
# R_seurat5 /net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/singleCell/scripts/R/seurat/integration/seurat_integration.R \
# "analysis/obj/individual/S283_Expression1.SeuratObj.QC.rds|analysis/obj/individual/S284_Expression1.SeuratObj.QC.rds|analysis/obj/individual/S285_Expression1.SeuratObj.QC.rds|analysis/obj/individual/S286_Expression1.SeuratObj.QC.rds|analysis/obj/individual/S287_Expression1.SeuratObj.QC.rds" \
# "analysis/obj/obj_seurat.SCT.5samples.rds" \
# "sct_CCA"


options(future.globals.maxSize = 150*1024 * 1024^2) #150GB
options(scipen = 999)

options(warn=-1)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(ggplot2)
})
options(warn=0)


# functions  ------------------------------------------------------------------------------------- ----
get_barcode_df=function(object_in){
  data.frame(barcode=row.names(object_in@meta.data),stringsAsFactors = F)
}



# temp parameters ------------------------------------------------------------------------------------- ----
# setwd("/net/nfs-irwrsrchnas01/labs/zgu_grp/Group/Grp_AncaPasca/omics/")
# 
files_rds="analysis/obj/individual/S283_Expression1.SeuratObj.QC.rds,
analysis/obj/individual/S284_Expression1.SeuratObj.QC.rds,
analysis/obj/individual/S285_Expression1.SeuratObj.QC.rds,
analysis/obj/individual/S286_Expression1.SeuratObj.QC.rds,
analysis/obj/individual/S287_Expression1.SeuratObj.QC.rds"
# 
# file_rds_out="analysis/obj/obj_seurat.SCT.5samples.rds"
# 
# method="sct_cca"

# parameters ------------------------------------------------------------------------------------- ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

files_rds=args[1]
file_rds_out=args[2]
method=args[3]

# method="sct_cca"
# method="harmony"

# other parameters ------------------------------------------------------------------------------------- ----
do_subset=F;n_max=2000

resolution_FindClusters=2
dims_used=50

w_umap=12;h_umap=10

# Anchor-based CCA integration (method=CCAIntegration)
# Anchor-based RPCA integration (method=RPCAIntegration)
# Harmony (method=HarmonyIntegration)
# FastMNN (method= FastMNNIntegration)
# scVI (method=scVIIntegration)

# pharse parameters ------------------------------------------------------------------------------------- ----
list_rds=unlist(strsplit(files_rds,split = "[,\n]"))
list_rds=list_rds[!list_rds==""]

dir_out=paste0(dirname(file_rds_out),"/")
if(!dir.exists(dir_out)){dir.create(dir_out,recursive = T)}

label=word(basename(file_rds_out),-2,sep="[.]")

prefix=gsub(".rds","",file_rds_out)

# read obj  ------------------------------------------------------------------------------------- ----
cat("reading obj...\n")

list_obj = lapply(list_rds, function(file_rds){
  obj_one=readRDS(file_rds)
  CreateSeuratObject(counts = GetAssayData(obj_one,assay = "RNA",layer = "count"),meta.data = obj_one@meta.data)
})


n_obj=length(list_obj)

cat(paste0(n_obj," objects in the integration list.\n"))


cat("Done Reading Obj\n")
if(do_subset){
  cat("subseting\n")
  for(i in 1:n_obj){
    set.seed(10)
    barcode_sample=get_barcode_df(list_obj[[i]]) 
    
    if(nrow(barcode_sample)>n_max){
      barcode_sample=barcode_sample %>% sample_n(n_max)
    } else {
      barcode_sample=barcode_sample
    }
    
    list_obj[[i]]=subset(list_obj[[i]],cells=barcode_sample$barcode)
  }
  cat("Done subseting\n")
}

# integration  ------------------------------------------------------------------------------------- ----
cat("running integration...\n")
obj_all <- merge(list_obj[[1]], 
                 y = list_obj[2:n_obj],
                 add.cell.ids = sapply(list_obj[1:n_obj], function(x) unique(x$group)), 
                 project = label)

table(obj_all$group)

# batch correctiong  ------------------------------------------------------------------------------------- ----

if (tolower(method)=="sct_cca"){
  cat("running batch integration (method=CCAIntegration, normalization.method = 'SCT')...\n")
  obj_all <- SCTransform(obj_all,verbose=F)
  obj_all <- RunPCA(obj_all,verbose = F)
  obj_all <- IntegrateLayers(object = obj_all, method = CCAIntegration, normalization.method = "SCT", verbose = F)
  
  obj_all <- RunUMAP(obj_all, dims = 1:dims_used, reduction = "integrated.dr")
  obj_all <- FindNeighbors(obj_all, reduction = "integrated.dr", dims = 1:dims_used)
  obj_all <- FindClusters(obj_all, resolution = resolution_FindClusters)
}


if (tolower(method)=="harmony"){
  cat("running batch integration (method=HarmonyIntegration, orig.reduction = 'pca')...\n")
  
  obj_all <- NormalizeData(obj_all,verbose = F)
  obj_all <- FindVariableFeatures(obj_all,verbose = F)
  obj_all <- ScaleData(obj_all,verbose = F)
  obj_all <- RunPCA(obj_all,verbose = F)
  
  obj_all <- IntegrateLayers(
    object = obj_all, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    verbose = FALSE
  )
  
  obj_all <- RunUMAP(obj_all, reduction = "harmony", dims = 1:dims_used, reduction.name = "umap.cca")
  
  obj_all <- FindNeighbors(obj_all, reduction = "harmony", dims = 1:dims_used)
  obj_all <- FindClusters(obj_all, resolution = resolution_FindClusters)
}


cat("Done integrating...\n")

# draw plot -------------------
DimPlot(obj_all,group.by = "group")
ggsave(paste0(prefix,"dimPlot.group.png"),width=w_umap,height = h_umap)
ggsave(paste0(prefix,"dimPlot.group.pdf"),width=w_umap,height = h_umap)

DimPlot(obj_all,group.by = "seurat_clusters",label = T,repel = T)
ggsave(paste0(prefix,"dimPlot.seurat_clusters.png"),width=w_umap,height = h_umap)
ggsave(paste0(prefix,"dimPlot.seurat_clusters.pdf"),width=w_umap,height = h_umap)

# save output ---------------------
cat("Saving output...\n")
saveRDS(obj_all,file=file_rds_out)
obj_all
cat("Done\n")





























