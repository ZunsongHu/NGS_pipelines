# source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

#temp parameters ---- ---- ---- ---- ---- ----
setwd("C:/Users/zuhu/Desktop/project/yuhongWang/cmv/")
# 
# rds="analysis/obj/CD34.SeuratObj.QC.rds"
dir_out="analysis/featurePlot/celltypist_features/"
# file_celltypeMarkers="C:/Users/zuhu/Desktop/script/R/singleCell/celltypist.celltypeMarkers.tsv"

#parameters ---- ---- ---- ---- ---- ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds=args[1]
dir_out=args[2]
file_celltypeMarkers=args[3]

#create dir_out -----
if(!dir.exists(dirname(dir_out))){dir.create(dirname(dir_out))}
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
obj_=readRDS(rds)
# obj_=obj_logC

gene_name_all=row.names(x=obj_)
# head(gene_name_all)

#read markers ----
df_celltypeMarkers=read_tsv(file_celltypeMarkers)

celltype=df_celltypeMarkers$celltype[1]

for (celltype in df_celltypeMarkers$celltype){
  markers=trimws(unlist(strsplit(df_celltypeMarkers$markers[df_celltypeMarkers$celltype==celltype],split = ",")))
  
  if(any(markers %in% gene_name_all)){
    i=sum(markers %in% gene_name_all)
    width=i*7
    height=6

    FeaturePlot_scCustom(seurat_object = obj_, features = markers, colors_use = viridis_plasma_dark_high, na_color = "lightgray")
    
    ggsave(paste0(dir_out,celltype,".png"),width=width,height=height)
    ggsave(paste0(dir_out,celltype,".pdf"),width=width,height=height)
    
  }
}

cat("Done\n")

