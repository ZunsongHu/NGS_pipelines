# set ---------------------------
# library(GenomicRanges)
# library(csaw)
source("/home/zgu_labs/bin/R/functions.R")
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

# testing parameters ----
# setwd("/scratch/zuhu/project/Grp_NadiaCarlesso/atac/")
# 
# dir_out=        "test/"
# label=          "test"
# file_samplelist=      "23.list"
# obj="analysis/peak_counts/obj.All.filtered.rds"


# parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

dir_out=          args[1]
label=            args[2]
file_samplelist=  args[3]
obj=              args[4]

# create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

# read sample list ----
samplelist=readLines(file_samplelist)

# read obj ----
obj_=readRDS(obj)

obj_out=obj_
obj_out$SE=obj_out$SE[,samplelist]

message(paste0(nrow(obj_out$SE)," peaks, ",ncol(obj_out$SE)," samples left after filtering"))

saveRDS(obj_out,paste0(dir_out,"obj.",label,".filtered",".rds"))


