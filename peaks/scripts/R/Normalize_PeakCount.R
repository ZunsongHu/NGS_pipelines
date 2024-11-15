library(dplyr)
library(stringr)
library(DESeq2)

gc(rm(list=ls()))

source("/home/zgu_labs/bin/R/functions.R")
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

#temp parameters ----
# setwd("/scratch/zuhu/project/ZhaohuiGu/PAX5/CUTTAG/")
# count_bed="/scratch/zuhu/project/ZhaohuiGu/PAX5/CUTTAG/out_raw/C13_C1/CutTag/macs2PeakMergeAll_readsCount/C13_C1.count_PeakMergeAll.bed,
# /scratch/zuhu/project/ZhaohuiGu/PAX5/CUTTAG/out_raw/C14_C1/CutTag/macs2PeakMergeAll_readsCount/C14_C1.count_PeakMergeAll.bed,
# /scratch/zuhu/project/ZhaohuiGu/PAX5/CUTTAG/out_raw/C15_C1/CutTag/macs2PeakMergeAll_readsCount/C15_C1.count_PeakMergeAll.bed,
# /scratch/zuhu/project/ZhaohuiGu/PAX5/CUTTAG/out_raw/C16_C1/CutTag/macs2PeakMergeAll_readsCount/C16_C1.count_PeakMergeAll.bed,
# /scratch/zuhu/project/ZhaohuiGu/PAX5/CUTTAG/out_raw/C17_C1/CutTag/macs2PeakMergeAll_readsCount/C17_C1.count_PeakMergeAll.bed,
# /scratch/zuhu/project/ZhaohuiGu/PAX5/CUTTAG/out_raw/C18_C1/CutTag/macs2PeakMergeAll_readsCount/C18_C1.count_PeakMergeAll.bed"
# file_anno="out_raw/stats/macs2PeakMergeAll/macs2PeakMergeAll.anno.tsv"
# control_id="C13_C1"
# dir_out="out_raw/stats/macs2PeakMergeAll/"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)
count_bed=args[1]
file_anno=args[2]
control_id=args[3]
dir_out=args[4]


#create output foler ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read count ----
count_bed_list=gsub("\n","",unlist(strsplit(count_bed,"[,]")))

for(count_bed_one in count_bed_list){
  COH_sample=word(basename(count_bed_one),1,sep="[.]")
  df_bedCount=read.table(count_bed_one,stringsAsFactors = F,sep="\t") %>%
    transmute(
      bedPos=paste0(V1,":",V2,"-",V3),
      length=V3-V2,
      count=V4,
    )
  names(df_bedCount)[3]=COH_sample
  
  if(match(count_bed_one,count_bed_list)==1){df_bedCount_all=df_bedCount}
  if(match(count_bed_one,count_bed_list)>=2){df_bedCount_all=df_bedCount_all %>% left_join(df_bedCount)}
}


#get vst ----
matrix_count=df_bedCount_all[3:ncol(df_bedCount_all)]
row.names(matrix_count)=df_bedCount_all$bedPos

df_vst=round(run_vst_matrix_df(matrix_count),4)
names(df_vst)=paste0(names(df_vst),".vst")

#get CPM ----
df_matrix_CPM=as.data.frame(round(DGEobj.utils::convertCounts(countsMatrix = as.matrix(matrix_count),unit = "CPM")),4)
names(df_matrix_CPM)=paste0(names(df_matrix_CPM),".CPM")
#get TPM ----
df_matrix_tpm=as.data.frame(round(DGEobj.utils::convertCounts(countsMatrix = as.matrix(matrix_count),unit = "TPM",geneLength=df_bedCount_all$length)),4)
names(df_matrix_tpm)=paste0(names(df_matrix_tpm),".TPM")

names(matrix_count)=paste0(names(matrix_count),".count")


#get comparision ----
control_vst_value=unlist(df_vst[paste0(control_id,".vst")])
df_matrix_vst_diff=round(df_vst-control_vst_value,4)
names(df_matrix_vst_diff)=paste0(names(df_matrix_vst_diff),".diff")

df_matrix_vst_foldChange=round(df_vst/control_vst_value,4)
names(df_matrix_vst_foldChange)=paste0(names(df_matrix_vst_foldChange),".foldChange")

#get df out ---
df_anno=read_tsv(file_anno) %>% select(bedPos,mutClass,dist2gene,transcript,gene,transStart,transEnd,exonPos,exonLen)

df_out=df_anno %>% left_join(bind_cols(df_bedCount_all[1:2],matrix_count,df_vst,df_matrix_CPM,df_matrix_tpm,df_matrix_vst_diff,df_matrix_vst_foldChange))

write_tsv(df_out,paste0(dir_out,"macs2PeakMergeAll_noramlization.tsv"))




