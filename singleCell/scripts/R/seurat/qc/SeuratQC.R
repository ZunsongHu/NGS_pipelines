# rm(list=ls())
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")
library(DoubletFinder)
source("/home/zgu_labs/bin/R/SingleCell/Seurat/doubleletFinder_hu.R")

# #temp parameters ----
# setwd("/scratch/zuhu/project/srivi/singlecell/")
# rds_in="out_raw/health_5p5K/Seurat/health_5p5K.SeuratObj.raw.rds"
# file_MTgeneList="NULL"
# dir_raw="test/"
# dir_qc="test/"
# rds_out="test/test.rds"
# species="mouse"
# transform="SCT"

# species="human"
# transform="ordinary"

#get parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=args[1]
file_MTgeneList=args[2]
dir_raw=args[3]

dir_qc=args[4]
rds_out=args[5]

species=args[6]
transform=args[7]

run_cellCycle=T

cutoff_MT=10

#read rds
obj_raw=readRDS(rds_in)

if(!dir.exists(dir_raw)){dir.create(dir_raw)}
if(!dir.exists(dir_qc)){dir.create(dir_qc)}

# #subset object
# df_barcode=get_barcode_df(obj_raw)
# set.seed(10)
# df_barcode_r=df_barcode %>% sample_n(1500)
# obj_raw=subset(obj_raw,cells=df_barcode_r$barcode)

if(species=="monkey"){species="human"}

#add QC measurements 
if(file_MTgeneList=="NULL"){
  if(species=="human"){
    obj_raw=add_qc_measurements(obj_raw,getMTgeneAuto = T)
  }
  
  if(species=="mouse"){
    obj_raw=add_qc_measurements(obj_raw,MTgene = T,housekeepinngGene = hkgenes_mouse)
  }
}

if(!file_MTgeneList=="NULL"){
  MTgeneList=readLines(file_MTgeneList)
  obj_raw=add_qc_measurements(obj_raw,getMTgeneAuto = F,MTgene = MTgeneList,housekeepinngGene = hkgenes_human)
}

#Normalization 
if(transform=="SCT"){
  #SC transform
  obj_raw=SCTransform(obj_raw)
}

if(transform=="ordinary"){
  #SC transform
  obj_raw=ScaleData(obj_raw)
  obj_raw=NormalizeData(obj_raw)
  obj_raw=FindVariableFeatures(obj_raw)
}

#run pca to cluster
obj_raw <- RunPCA(obj_raw, verbose = FALSE)
obj_raw <- RunUMAP(obj_raw, dims = 1:50, verbose = FALSE)
obj_raw <- FindNeighbors(obj_raw, dims = 1:50, verbose = FALSE)
obj_raw <- FindClusters(obj_raw, verbose = FALSE,resolution = 1)

#run doublelets
cat("####################################################################################################")
cat("remove doublets")
cat("####################################################################################################")

obj_raw=get_doublet_DoubletFinder(obj_raw)

DimPlot(obj_raw,group.by=names(obj_raw@meta.data)[grepl("DF.classifications",names(obj_raw@meta.data))])+
  FeaturePlot(obj_raw,features = names(obj_raw@meta.data)[grepl("pANN",names(obj_raw@meta.data))])+plot_layout(ncol=2)

ggsave(paste0(dir_raw,"QC_DoubletFinder_raw",".png"),width=12,height = 7)
ggsave(paste0(dir_raw,"QC_DoubletFinder_raw",".pdf"),width=12,height = 7)

#draw plot
FeaturePlot(obj_raw,features=c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribosomal"),ncol = 2)
ggsave(paste0(dir_raw,"QC_basic_raw_FeaturePlot",".png"),width=12,height = 12)
ggsave(paste0(dir_raw,"QC_basic_raw_FeaturePlot",".pdf"),width=12,height = 12)

VlnPlot(obj_raw,features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribosomal"),ncol = 2,group.by = "group")
ggsave(paste0(dir_raw,"QC_basic_raw_VlnPlot",".png"),width=12,height = 12)
ggsave(paste0(dir_raw,"QC_basic_raw_VlnPlot",".pdf"),width=12,height = 12)

FeatureScatter(obj_raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "group")+
  FeatureScatter(obj_raw, feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "group")+
  FeatureScatter(obj_raw, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "group")

ggsave(paste0(dir_raw,"QC_basic_raw_FeatureScatter",".png"),width=18,height = 6)
ggsave(paste0(dir_raw,"QC_basic_raw_FeatureScatter",".pdf"),width=18,height = 6)


if(run_cellCycle==T){
  #get cell cycle
  obj_raw = CellCycleScoring(object = obj_raw, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
  
  FeaturePlot(obj_raw,features="S.Score")+
    FeaturePlot(obj_raw,features="G2M.Score")+
    DimPlot(obj_raw,group.by='Phase')+plot_layout(ncol=2)
  
  ggsave(paste0(dir_raw,"QC_CellCycle_raw",".png"),width=12,height = 12)
  ggsave(paste0(dir_raw,"QC_CellCycle_raw",".pdf"),width=12,height = 12)
}

df_singlelet=get_metadataVar(obj_raw,names(obj_raw@meta.data)[grepl("DF.classifications",names(obj_raw@meta.data))])
df_singlelet=df_singlelet[df_singlelet[,2]=="Singlet",]
write_tsv(df_singlelet,paste0(dir_raw,"df_singlelet.tsv"))

#Remove doublets----------------------------------------------------------------
obj_afterRemovingDoublet=subset(obj_raw,cells=df_singlelet$barcode)

cat("####################################################################################################")
cat(paste0("number of cells removed because of doublelets: ",ncol(obj_raw)-ncol(obj_afterRemovingDoublet)))
cat("####################################################################################################")

#add QC measurements 
if(file_MTgeneList=="NULL"){
  if(species=="human"){
    obj_afterRemovingDoublet=add_qc_measurements(obj_afterRemovingDoublet)
  }
  if(species=="mouse"){
    obj_afterRemovingDoublet=add_qc_measurements(obj_afterRemovingDoublet,MTgene = T,housekeepinngGene = hkgenes_mouse)
  }
}

if(!file_MTgeneList=="NULL"){
  MTgeneList=readLines(file_MTgeneList)
  obj_afterRemovingDoublet=add_qc_measurements(obj_afterRemovingDoublet,getMTgeneAuto = F,
                                               MTgene = MTgeneList,housekeepinngGene = hkgenes_human)
}

if(transform=="SCT"){
  #SC transform
  obj_afterRemovingDoublet=SCTransform(obj_afterRemovingDoublet)
}

if(transform=="ordinary"){
  #SC transform
  obj_afterRemovingDoublet=NormalizeData(obj_afterRemovingDoublet)
  obj_afterRemovingDoublet=FindVariableFeatures(obj_afterRemovingDoublet)
  obj_afterRemovingDoublet=ScaleData(obj_afterRemovingDoublet)
}

obj_afterRemovingDoublet <- RunPCA(obj_afterRemovingDoublet, verbose = FALSE)
obj_afterRemovingDoublet <- RunUMAP(obj_afterRemovingDoublet, dims = 1:50, verbose = FALSE)
obj_afterRemovingDoublet <- FindNeighbors(obj_afterRemovingDoublet, dims = 1:50, verbose = FALSE)
obj_afterRemovingDoublet <- FindClusters(obj_afterRemovingDoublet, resolution = 1,verbose = FALSE)

# obj_afterRemovingDoublet=get_doublet_DoubletFinder(obj_afterRemovingDoublet)

FeaturePlot(obj_afterRemovingDoublet,features=c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribosomal"),ncol = 2)
ggsave(paste0(dir_raw,"QC_basic_afterRemovingDoublet_FeaturePlot",".png"),width=12,height = 12)
ggsave(paste0(dir_raw,"QC_basic_afterRemovingDoublet_FeaturePlot",".pdf"),width=12,height = 12)

VlnPlot(obj_afterRemovingDoublet,features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribosomal"),ncol = 2,group.by = "group")
ggsave(paste0(dir_raw,"QC_basic_afterRemovingDoublet_VlnPlot",".png"),width=12,height = 12)
ggsave(paste0(dir_raw,"QC_basic_afterRemovingDoublet_VlnPlot",".pdf"),width=12,height = 12)


FeatureScatter(obj_afterRemovingDoublet, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "group")+
  FeatureScatter(obj_afterRemovingDoublet, feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "group")+
  FeatureScatter(obj_afterRemovingDoublet, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "group")

ggsave(paste0(dir_raw,"QC_basic_afterRemovingDoublet_FeatureScatter",".png"),width=18,height = 6)
ggsave(paste0(dir_raw,"QC_basic_afterRemovingDoublet_FeatureScatter",".pdf"),width=18,height = 6)

if(run_cellCycle==T){
  obj_afterRemovingDoublet = CellCycleScoring(object = obj_afterRemovingDoublet, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
  
  FeaturePlot(obj_afterRemovingDoublet,features="S.Score")+
    FeaturePlot(obj_afterRemovingDoublet,features="G2M.Score")+
    DimPlot(obj_afterRemovingDoublet,group.by='Phase')+plot_layout(ncol=2)
  
  ggsave(paste0(dir_raw,"QC_CellCycle_afterRemovingDoublet_FeaturePlot",".png"),width=12,height = 12)
  ggsave(paste0(dir_raw,"QC_CellCycle_afterRemovingDoublet_FeaturePlot",".pdf"),width=12,height = 12)
}

#Remove low quality cells ----------------------------------------------------------------
df_cutoff=get_cutoff(obj_afterRemovingDoublet)
write_tsv(df_cutoff,paste0(dir_qc,"df_cutoff.tsv"))

obj_QC1=subset(obj_afterRemovingDoublet, subset = percent.mt < cutoff_MT)
cat("####################################################################################################")
cat(paste0("number of cells removed because of dead cells: ",ncol(obj_afterRemovingDoublet)-ncol(obj_QC1)))
cat("####################################################################################################")

obj_QC=subset(obj_QC1, subset = nFeature_RNA > 200 & nFeature_RNA < df_cutoff$medianAdd3MAD[df_cutoff$var=="nFeature_RNA"])
cat("####################################################################################################")
cat(paste0("number of cells removed because of outlier cells of nFeature: ",ncol(obj_QC1)-ncol(obj_QC)))
cat("####################################################################################################")

#add QC measurements
if(file_MTgeneList=="NULL"){
  if(species=="human"){
    obj_QC=add_qc_measurements(obj_QC)
  }
  
  if(species=="mouse"){
    obj_QC=add_qc_measurements(obj_QC,MTgene = T,housekeepinngGene = hkgenes_mouse)
  }
}
  
if(!file_MTgeneList=="NULL"){
  MTgeneList=readLines(file_MTgeneList)
  obj_QC=add_qc_measurements(obj_QC,getMTgeneAuto = F,
                             MTgene = MTgeneList,housekeepinngGene = hkgenes_human)
}
  
if(transform=="SCT"){
  #SC transform
  obj_QC=SCTransform(obj_QC)
}

if(transform=="ordinary"){
  #SC transform
  obj_QC=NormalizeData(obj_QC)
  obj_QC=FindVariableFeatures(obj_QC)
  obj_QC=ScaleData(obj_QC)
}

dims=50
n.neighbors=5
obj_QC <- RunPCA(obj_QC, verbose = FALSE,npcs = dims)
obj_QC <- RunUMAP(obj_QC, dims = 1:dims, verbose = FALSE,n.neighbors = n.neighbors)
obj_QC <- FindNeighbors(obj_QC, dims = 1:dims, verbose = FALSE)
obj_QC <- FindClusters(obj_QC, resolution = 1,verbose = FALSE)


FeaturePlot(obj_QC,features=c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribosomal"),ncol = 2)
ggsave(paste0(dir_qc,"QC_basic_afterRemovingDoublet_remoingLowQuatlity_FeaturePlot",".png"),width=12,height = 12)
ggsave(paste0(dir_qc,"QC_basic_afterRemovingDoublet_remoingLowQuatlity_FeaturePlot",".pdf"),width=12,height = 12)

VlnPlot(obj_QC,features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribosomal"),ncol = 2,group.by = "group")
ggsave(paste0(dir_qc,"QC_basic_afterRemovingDoublet_remoingLowQuatlity_VlnPlot",".png"),width=12,height = 12)
ggsave(paste0(dir_qc,"QC_basic_afterRemovingDoublet_remoingLowQuatlity_VlnPlot",".pdf"),width=12,height = 12)


FeatureScatter(obj_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",group.by = "group")+
  FeatureScatter(obj_QC, feature1 = "nFeature_RNA", feature2 = "percent.mt",group.by = "group")+
  FeatureScatter(obj_QC, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "group")

ggsave(paste0(dir_qc,"QC_basic_afterRemovingDoublet_remoingLowQuatlity_FeatureScatter",".png"),width=18,height = 6)
ggsave(paste0(dir_qc,"QC_basic_afterRemovingDoublet_remoingLowQuatlity_FeatureScatter",".pdf"),width=18,height = 6)

if(run_cellCycle==T){
  
  FeaturePlot(obj_QC,features="S.Score")+
    FeaturePlot(obj_QC,features="G2M.Score")+
    DimPlot(obj_QC,group.by='Phase')+plot_layout(ncol=2)
  
  ggsave(paste0(dir_qc,"QC_CellCycle_afterRemovingDoublet_remoingLowQuatlity_FeaturePlot",".png"),width=12,height = 12)
  ggsave(paste0(dir_qc,"QC_CellCycle_afterRemovingDoublet_remoingLowQuatlity_FeaturePlot",".pdf"),width=12,height = 12)
  obj_QC = CellCycleScoring(object = obj_QC, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
}


df_bc=get_barcode_df(obj_QC)
write.table(df_bc,paste0(gsub(".rds","",rds_out),".barcode.list"),row.names = F,col.names = F,quote = F)

saveRDS(obj_QC,file=rds_out)



