library(Seurat)

gc(rm(list=ls()))
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")


# parameters temp ----
# setwd("/scratch/zuhu/tmp/SC_CMV/")
# scid="NY"
# rds_in="out_raw/NY/Seurat/NY.SeuratObj.raw.rds"
# rds_out="out_raw/NY/Seurat/NY.SeuratObj.QC.rds"
# dir_out="out_raw/NY/Seurat/QCplot/"
# transform="SCT"

run_cellCycle=T

# parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

scid=args[1]
rds_in=args[2]
rds_out=args[3]
dir_out=args[4]
transform=args[5]

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read rds ----
obj_=readRDS(rds_in)

# obj_=subset(obj_,cells=(get_barcode_df(obj_) %>% sample_n(1000))$barcode)

obj_[["group"]]=scid

#run normalization ----
if(transform=="SCT"){
  #SC transform
  obj_=SCTransform(obj_)
}

if(transform=="ordinary"){
  #SC transform
  obj_=ScaleData(obj_)
  obj_=NormalizeData(obj_)
  obj_=FindVariableFeatures(obj_)
}

#run pca to cluster ----
obj_ <- RunPCA(obj_, verbose = FALSE)
obj_ <- RunUMAP(obj_, dims = 1:50, verbose = FALSE)
obj_ <- FindNeighbors(obj_, dims = 1:50, verbose = FALSE)
obj_ <- FindClusters(obj_, verbose = FALSE)


#draw plot ----
DimPlot(obj_,group.by = "seurat_clusters",label=T,repel = T)
ggsave(paste0(dir_out,"DimPlot_SeuratClusters",".png"),width=8,height = 8)
ggsave(paste0(dir_out,"DimPlot_SeuratClusters",".pdf"),width=8,height = 8)

FeaturePlot(obj_,features=c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribosomal"),ncol = 2)
ggsave(paste0(dir_out,"FeaturePlot_QCbasic",".png"),width=12,height = 12)
ggsave(paste0(dir_out,"FeaturePlot_QCbasic",".pdf"),width=12,height = 12)

VlnPlot(obj_,features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribosomal"),ncol = 2,group.by = "group")
ggsave(paste0(dir_out,"VlnPlot_QCbasic",".png"),width=12,height = 12)
ggsave(paste0(dir_out,"VlnPlot_QCbasic",".pdf"),width=12,height = 12)

FeatureScatter(obj_, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "group")+
  FeatureScatter(obj_, feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "group")+
  FeatureScatter(obj_, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "group")
ggsave(paste0(dir_out,"FeatureScatter_QCbasic",".png"),width=18,height = 6)
ggsave(paste0(dir_out,"FeatureScatter_QCbasic",".pdf"),width=18,height = 6)

if(run_cellCycle==T){
  #get cell cycle
  obj_ = CellCycleScoring(object = obj_, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
  
  FeaturePlot(obj_,features="S.Score")+
    FeaturePlot(obj_,features="G2M.Score")+
    DimPlot(obj_,group.by='Phase')+plot_layout(ncol=2)
  
  ggsave(paste0(dir_out,"FeaturePlot_CellCycle",".png"),width=12,height = 12)
  ggsave(paste0(dir_out,"FeaturePlot_CellCycle",".pdf"),width=12,height = 12)
}

saveRDS(obj_,file=rds_out)





























































