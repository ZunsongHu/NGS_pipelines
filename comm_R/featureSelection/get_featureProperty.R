#rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
# library(SummarizedExperiment, quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# obj_rds="../tsne/2042_middleGroup/obj_2042MiddleGroup_totalRNA.rds"
# label=          "2042_middleGroup_totalRNA"
# dir_out=        "test/"


#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=        args[1]
label=          args[2]
dir_out=        args[3]


high_trim_n=1
low_trim_n=1
distance_abs=1.5

#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}


#read obj ----
cat(paste0("Reading obj ....","\n"))
obj_=readRDS(obj_rds)
df_vst_all=assays(obj_$SE)[["vst"]] 

# df_vst_all=df_vst_all %>% sample_n(1000)

#get median values ----
df_meidanDiag=df_vst_all %>% mutate(feature_id=row.names(df_vst_all)) %>% 
  reshape2::melt() %>% 
  mutate(COH_sample=variable) %>% 
  left_join(as.data.frame(colData(obj_$SE))) %>% 
  group_by(feature_id,diag) %>% 
  mutate(median_value=median(value)) %>% 
  select(feature_id,diag,median_value) %>% 
  distinct() %>% ungroup() %>% 
  spread(diag,median_value)



#get label for outlier values ----
diag_levles=names(df_meidanDiag)[-1]

get_outlierLabel=function(x){
  n_=length(x)
  
  df_x=data.frame(value=x,name=names(x),stringsAsFactors = F) 
  
  df_x1=df_x %>% 
    arrange(value,name) %>% slice_head(n=(n_-high_trim_n)) %>% 
    arrange(desc(value),name) %>% slice_head(n=(n_-(high_trim_n+low_trim_n)))
  
  mean_=mean(df_x1$value)
  sd_=sd(df_x1$value)
  
  df_x2=df_x %>% 
    mutate(
      meanMinus1SD=mean_-sd_,meanPlus1SD=mean_+sd_,
      meanMinus2SD=mean_-2*sd_,meanPlus2SD=mean_+2*sd_,
      meanMinus3SD=mean_-3*sd_,meanPlus3SD=mean_+3*sd_,
      # meanMinus4SD=mean_-4*sd_,meanPlus4SD=mean_+4*sd_,
      # meanMinus6SD=mean_-6*sd_,meanPlus6SD=mean_+6*sd_,
      # meanMinus10SD=mean_-10*sd_,meanPlus10SD=mean_+10*sd_,
      foldChange2Mean=round(value/mean_,4),
      distance2Mean=round(value-mean_,4),
      
      mean1SD_picked=ifelse((value >=meanPlus1SD | value <= meanMinus1SD) ,name,NA),
      mean2SD_picked=ifelse((value >=meanPlus2SD | value <= meanMinus2SD) ,name,NA),
      mean3SD_picked=ifelse((value >=meanPlus3SD | value <= meanMinus3SD) ,name,NA)
    ) %>% 
    filter(abs(distance2Mean) >= distance_abs)
  
  paste0(
    paste0(paste0(paste0(df_x2$mean1SD_picked[!is.na(df_x2$mean1SD_picked)],"(",df_x2$foldChange2Mean[!is.na(df_x2$mean1SD_picked)],",",df_x2$distance2Mean[!is.na(df_x2$mean1SD_picked)],")"),collapse = ";")),"|",
    paste0(paste0(paste0(df_x2$mean2SD_picked[!is.na(df_x2$mean2SD_picked)],"(",df_x2$foldChange2Mean[!is.na(df_x2$mean2SD_picked)],",",df_x2$distance2Mean[!is.na(df_x2$mean2SD_picked)],")"),collapse = ";")),"|",
    paste0(paste0(paste0(df_x2$mean3SD_picked[!is.na(df_x2$mean3SD_picked)],"(",df_x2$foldChange2Mean[!is.na(df_x2$mean3SD_picked)],",",df_x2$distance2Mean[!is.na(df_x2$mean3SD_picked)],")"),collapse = ";"))
  )
}

get_foldchage_=function(x){ifelse(!grepl(";",x),as.numeric(word(word(x,2,sep="[(]"),1,sep=",")),NA)}
get_distance_=function(x){ifelse(!grepl(";",x),gsub("[)]","",word(word(x,2,sep="[(]"),2,sep=",")),NA)}

# 
# 
# x=apply(df_meidanDiag[1:1000,diag_levles],1,get_outlierLabel)
# 
# df_=data.frame(label_outlier=x,stringsAsFactors = F) %>%
#   separate(label_outlier,sep="[|]",into=c(c("mean1SD_picked","mean2SD_picked","mean3SD_picked"))) %>% 
#   mutate(
#     mean1SD_picked=ifelse(mean1SD_picked=="(,)","",mean1SD_picked),
#     mean2SD_picked=ifelse(mean2SD_picked=="(,)","",mean2SD_picked),
#     mean3SD_picked=ifelse(mean3SD_picked=="(,)","",mean3SD_picked),
#     mean3SDPicked_foldchange2Mean=get_foldchage_(mean3SD_picked),
#     mean3SDPicked_distance2Mean=get_distance_(mean3SD_picked),
#     mean3SDPicked_distance2Mean=as.numeric(mean3SDPicked_distance2Mean)
#   ) 
  

df_meidanDiag$label_outlier=apply(df_meidanDiag[diag_levles],1,get_outlierLabel)

df_meidanDiag1=df_meidanDiag %>% 
  separate(label_outlier,sep="[|]",into=c(c("mean1SD_picked","mean2SD_picked","mean3SD_picked"))) %>% 
  mutate(
    mean1SD_picked=ifelse(mean1SD_picked=="(,)","",mean1SD_picked),
    mean2SD_picked=ifelse(mean2SD_picked=="(,)","",mean2SD_picked),
    mean3SD_picked=ifelse(mean3SD_picked=="(,)","",mean3SD_picked),
    mean3SD_picked_singleton=ifelse(!grepl(";",mean3SD_picked),word(mean3SD_picked,1,sep="[()]"),NA),
    mean3SD_picked_singleton=ifelse(mean3SD_picked_singleton=="",NA,mean3SD_picked_singleton),
    mean3SDPicked_foldchange2Mean=get_foldchage_(mean3SD_picked),
    mean3SDPicked_distance2Mean=get_distance_(mean3SD_picked),
    mean3SDPicked_distance2Mean=as.numeric(mean3SDPicked_distance2Mean)
  ) %>% 
  arrange(mean3SD_picked_singleton,mean3SDPicked_distance2Mean)


write_tsv(df_meidanDiag1,paste0(dir_out,"df_featureProperty.",label,".tsv"))





