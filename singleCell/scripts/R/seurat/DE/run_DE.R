source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

#temp parameters ----
# setwd("/scratch/zuhu/project/srivi/singlecell/")
# 
# file_obj="analysis/obj/obj.11samples.logNormalizationBatchCorrected.selected.rds"
# dir_out="test/"
# test_level="B cell"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_obj=args[1]
dir_out=args[2]
test_level=args[3]

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
cat("Reading obj...\n")
obj_=readRDS(file_obj)

#subset ----
set.seed(10)
obj_=subset(obj_, downsample = 200)

#run DE ----
cat("Running DE analysis...\n")
df_marker = FindMarkers(obj_, ident.1 = test_level, ident.2 = NULL, only.pos = TRUE)
df_marker=bind_cols(
  data.frame(gene_name=row.names(df_marker),stringsAsFactors = F),
  df_marker
)

write_tsv(df_marker,paste0(dir_out,"df_DE.",gsub("[/ ]","_",test_level),".tsv"))

cat("Done\n")


