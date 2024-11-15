# rm(list=ls())
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

# #temp parameters ----
setwd("/scratch/zuhu/project/srivi/singlecell/")
rds_in="out_raw/health_5p5K/Seurat/health_5p5K.SeuratObj.raw.rds"
rds_out="test/test.rds"
n=1500

#get parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=args[1]
rds_out=args[2]
n=as.numeric(args[3])

#create dir out ----
if(!dir.exists(dirname(rds_out))){dir.create(dirname(rds_out))}

#read rds
cat("reading obj ...\n")
obj_raw=readRDS(rds_in)

cat("Raw object:\n")
obj_raw
cat("\n")

#get coding gene names ----
df_bc=get_barcode_df(obj_raw)

set.seed(10)
df_bc=df_bc %>% sample_n(size = n)

obj_raw1=subset(obj_raw,cells=df_bc$barcode)

cat("New object:\n")
obj_raw1
cat("\n")
cat("Saving output object ...\n")
saveRDS(obj_raw1,rds_out)

cat("Done\n")

