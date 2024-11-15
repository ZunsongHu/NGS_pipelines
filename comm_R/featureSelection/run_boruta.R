options(future.globals.maxSize =10*1024 * 1024^2) #10GB

rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

library (SummarizedExperiment, quietly = T)
library(UBL, quietly = T)
library(Boruta)
library(stringr)

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=          args[1]
dataset_label=    args[2]
var_group=        args[3]
num.threads=      as.numeric(args[4])
dir_out=          args[5]

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# 
# obj_rds=          "2042_junctions/smote/1/smote_2042_junctions_diag_EachGroupN10.rds"
# dataset_label =   "2042_junctions"
# var_group=        "diag"
# dir_out=          "test/"
# num.threads=2



obj_type=         "matrix"
# obj_type=         "list"

features=NULL

#create dir out ----
if(!dir.exists(dirname(dir_out))){dir.create(dirname(dir_out))}
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
obj_=readRDS(obj_rds)

if(tolower(obj_type) %in% "matrix"){
  if(is.null(features)){
    features_in=names(obj_)[-match(var_group,names(obj_))]
    }
  if(!is.null(features)){
    # features_in=features
    features_in=names(obj_)[word(names(obj_),1,sep=":") %in% features]
  }
  
  value_group=unlist(obj_[var_group])
  df_in=obj_[features_in]
  df_in$group=value_group
}


if(tolower(obj_type) %in% "list"){
  if(is.null(features)){features_in=rownames(obj_$SE)}
  if(!is.null(features)){
    features_in=rownames(obj_$SE)[word(rownames(obj_$SE),1,sep=":") %in% features]}
  
  df_in=get_embeding_feature(obj_in = obj_,assay_name_in = "vst",features = c(var_group,features_in)) 
  names(df_in)[1]="group"
  df_in$group=as.factor(df_in$group)
}

dim(df_in)

#Run Boruta ----
print("Run Boruta")
set.seed(0)
boruta_train = Boruta(group ~ ., data = df_in, doTrace = 2,maxRuns=500,num.threads=num.threads)

#get output data frame -----
boruta_out_Tentative = TentativeRoughFix(boruta_train)

importance_boruta=attStats(boruta_train)
importance_boruta$gene_id_raw=row.names(importance_boruta)

importance_boruta_Tentative=attStats(boruta_out_Tentative)
names(importance_boruta_Tentative)=paste0(names(importance_boruta_Tentative),"_Tentative")
importance_boruta_Tentative$gene_id_raw=row.names(importance_boruta_Tentative)

df_boruta_out=importance_boruta %>% left_join(importance_boruta_Tentative) %>% 
  mutate(gene_id=gsub("`","",word(gene_id_raw,1,sep=":"))
)

df_out=info_gtf %>% 
  left_join(df_boruta_out) %>%
  filter(!is.na(meanImp))

#save output ----
obj_out=list(fit_boruta=boruta_train,
             boruta_out=df_boruta_out) 

print("Save to output")
saveRDS(obj_out,paste0(dir_out,"Boruta_",dataset_label,"_",var_group,".rds"))

# x=readRDS("2042_featureCountsExon/Boruta/2/Boruta_2042_featureCountsExon_EachGroupN10_diag.rds")
# 
# xx=x$boruta_out


