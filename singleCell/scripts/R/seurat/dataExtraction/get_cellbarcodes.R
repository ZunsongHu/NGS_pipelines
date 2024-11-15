#set ----------------------
gc(rm(list=ls()))
library(stringr)
library(reshape2)
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

setwd("/scratch/zuhu/project/ZhaohuiGu/SingleCellSeq/")

#temp parameters ----
rds_list="4.obj.list"
list_out="idlist/COHSC1_barcode.list"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_list=args[1]
list_out=args[2]

#read obj ----
files_list=readLines(rds_list)

file_one=files_list[1]

df_barcode=bind_rows(
  lapply(files_list, function(file_one){
    obj_1=readRDS(file_one)
    get_barcode_df(obj_1)
  })
)

write_tsv(df_barcode,list_out)






