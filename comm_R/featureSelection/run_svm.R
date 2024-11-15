options(future.globals.maxSize =10*1024 * 1024^2) #10GB

rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

library (SummarizedExperiment, quietly = T)
library(stringr)
library(caret)
library(kernlab)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# 
# obj_rds=          "/scratch/zuhu/project/BALL/tsne/1821/obj_1821.rds"
# dir_out=          "test/"
# dataset_label =   "1821"
# feature_panel=    "genes_1366_coding"
# featureN=         50
# obj_for_testing=  "../tsne/2955/obj_2955_forPrediction.rds"

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=          args[1]
dir_out=          args[2]
dataset_label=    args[3]
feature_panel=    args[4]
featureN=         as.numeric(args[5])
obj_for_testing=  args[6]
method=           args[7]

#create dir out ----
if(!dir.exists(dirname(dir_out))){dir.create(dirname(dir_out))}
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj for training ----
cat("Reading trainning object...\n")
obj_=readRDS(obj_rds)
# n_g=as.numeric(gsub("EachGroupN","",word(gsub(".rds","",basename(obj_rds)),-1,sep="_")))

#get feature N ----
featureN_total=length(obj_[[feature_panel]])

cat(paste0("Total feature N= ",featureN_total,"\n"))

featureN_top=ceiling(featureN_total/100)*100

if(featureN>featureN_top){cat("Feature N more than 100+total number of features in panel, no need to run\n")}

if(featureN<=featureN_top){
  featureN_=min(c(featureN_total,featureN))
  cat(paste0("Used feature N=",featureN_,"\n"))
  
  #read obj for testing ----
  cat("Reading test object...\n")
  obj_testing=readRDS(obj_for_testing)
  
  df_testing_all=as.data.frame(t(assays(obj_testing$SE)[["vst"]]))
  
  # str(df_testing_all)
  dim(df_testing_all)
  
  df_id=data.frame(COH_sample=row.names(df_testing_all),stringsAsFactors = F) %>% arrange(COH_sample)
  
  
  #get trainning matrix ----
  features_in=obj_[[feature_panel]][1:featureN_]
  
  df_in=get_features_df(obj_in = obj_,assay_name_in = "vst",features = c("diag",features_in)) 

  cat(paste0("n_features in matrix for training = ",dim(df_in)[2],"\n"))
  cat(paste0("n_samples in matrix for training = ",dim(df_in)[1],"\n"))
  
  #training svm ----
  df_train_=df_in[c("diag",features_in)]
  
  # df_train_$diag=ifelse( df_train_$diag=="Ph","Ph","Other")
  # table( df_train_$diag)
  # cat("Training model (unWeighted)...\n")
  # set.seed(10)
  # svm_fit = train(diag ~., data = df_train_,trControl = trainControl(method = "none"),preProcess = c("center","scale"),method = "svmLinear")
  

  cat("Training model...\n")
  # df_weights=df_train_["diag"] %>% group_by(diag) %>% mutate(n_g=n(),weight=round(100/n_g,4))
  # 
  set.seed(10)
  # ctrl = trainControl(method = "repeatedcv",number = 5,repeats = 5,sampling="smote")
  
  # ctrl = trainControl(method = "none") #none
  # ctrl = trainControl(method = "none", sampling="smote") #none + smote
  # 
  # ctrl = trainControl(method = "boot",number = 10) #boot
  # ctrl = trainControl(method = "boot",number = 10, sampling="smote") #boot + smote
  # 
  ctrl = trainControl(method = "repeatedcv",number = 10,repeats = 3,trim=TRUE,returnData = FALSE ) #repeatedcv
  # ctrl = trainControl(method = "repeatedcv",number = 10,repeats = 5, sampling="smote") #repeatedcv + smote

  set.seed(10)
  # svm_fit = train(diag ~., data = df_train_,preProcess = c("center","scale"),trControl=ctrl,method = "svmLinear3") #default one
  # method = "svmLinear3"
  svm_fit = train(diag ~., data = df_train_,preProcess = c("center","scale"),trControl=ctrl,method = method)
  
  
  cat("Predicting...\n")
  df_test_=df_testing_all[,features_in]
  pred_svm=caret::predict.train(svm_fit,df_test_)
    

  cat("output...\n")
  df_out=df_id %>% mutate(dataset_label=dataset_label,feature_panel=feature_panel,featureN=featureN_,pred_svm=pred_svm)

  #output ----
  file_out_base=paste0(dir_out,dataset_label,".",feature_panel,".FeatureN",featureN,".svm")
    
  cat(paste0(file_out_base,".pred.tsv\n"))
  write_tsv(df_out,paste0(file_out_base,".pred.tsv"))
    
  saveRDS(svm_fit,paste0(file_out_base,".model.rds"))
    
}





