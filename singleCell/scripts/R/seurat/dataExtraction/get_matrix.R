options(future.globals.maxSize = 10*1024 * 1024^2) #2GB

gc(rm(list=ls()))
library(Seurat)
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

# parameters temp ----
setwd("/scratch/zuhu/project/Grp_Swaminathan/singlecell/")

dir_matrix= "out_raw/A2_Expression1/Cellranger_Count_Expectnumber/filtered_feature_bc_matrix/"
id=         'A2_Expression1'
dir_out=    'out_raw/A2_Expression1/expMatrix/'

# parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

dir_matrix= args[1]
id=         args[2]
dir_out=    args[3]


#create output folder -------------------------
if(!dir.exists(dir_out)){dir.create(dir_out)}

#create seurat object 
message("Creating seurat object ...")
data = Read10X(data.dir = dir_matrix,gene.column=2)
obj_=CreateSeuratObject(counts = data)

#run normalization
message("Running normalization ...")
obj_=NormalizeData(obj_)

#get matrix
message("Getting matrix ...")
matrix_count=GetAssayData(obj_,assay = "RNA",slot = "counts")
matrix_normalization=round(as.matrix(GetAssayData(obj_,assay = "RNA",slot = "data")),4)

matrix_count=bind_cols(
  data.frame(feature=row.names(matrix_count),stringsAsFactors = F),
  as.data.frame(matrix_count)
)

matrix_normalization=bind_cols(
    data.frame(feature=row.names(matrix_normalization),stringsAsFactors = F),
    as.data.frame(matrix_normalization)
)

#Writing matrix
message("Writing matrix ...")
write_tsv(matrix_count,paste0(dir_out,id,".countMatrix.tsv"))
write_tsv(matrix_normalization,paste0(dir_out,id,".normalizationMatrix.tsv"))

message("Done getMatix")











