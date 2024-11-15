# options(future.globals.maxSize =10*1024 * 1024^2) #10GB
# gc(rm(list=ls()))
# rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

suppressWarnings(suppressMessages(require(SummarizedExperiment)))

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=            args[1]
dataset_label=      args[2]
var_group=          args[3]
file_feature_fixed= args[4]
seed_in=            as.numeric(args[5])
n_for_features=     as.numeric(args[6])
run_GA=             args[7]
dir_out=            args[8]

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# 
# obj_rds=            "2042/smote/smote_2042_diag_EachGroupN100.rds"
# dataset_label =     "2042_EachGroupN50"
# var_group=          "diag"
# file_feature_fixed= "2042/GA_summarise/GA_2042_EachGroupN100_NTrain50_NFeatures800_NKNN5_DimsPca50_top.list"
# seed_in=            1
# n_for_features=     800
# dir_out=            "test/"
# run_GA="Yes"



n_for_train=50
KNN_K=5
pca_dims=50
obj_type=         "matrix"
mutation_rate=0.01
features=gene_in

# features=features[1:100]
#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
obj_=readRDS(obj_rds)

if(tolower(obj_type)=="matrix"){
  value_group=unlist(obj_[var_group])
  df_in=obj_[c(features)]
  df_in=df_in[apply(df_in, 1, max) >= 5,]
  df_in$group=value_group
  df_in$obs=1:nrow(df_in)
}

#get train test dataset ----
n_for_features_=n_for_features-length(readLines(file_feature_fixed))
features_=features[!features %in% readLines(file_feature_fixed)]
set.seed(seed_in)
features_in=c(readLines(file_feature_fixed),sample(features,n_for_features_))

if(tolower(run_GA)=="yes"){
  n_fixed_in=length(readLines(file_feature_fixed))-ceiling(length(readLines(file_feature_fixed))*mutation_rate)
  set.seed(seed_in)
  features_fixed_in=sample(readLines(file_feature_fixed),n_fixed_in)
  
  features_=features[!features %in% features_fixed_in]
  n_for_features_=n_for_features-length(features_fixed_in)
  set.seed(seed_in)
  features_in=sort(unique(c(features_fixed_in,sample(features_,n_for_features_))))
}

set.seed(seed_in)
df_in=df_in[c("group","obs",features_in)]
df_train=df_in %>% group_by(group) %>% sample_n(n_for_train)
df_test=df_in[!df_in$obs %in% df_train$obs,]

#Run PCA ----
print("Run PCA")

pca_fit=prcomp(df_train[features_in], rank = pca_dims)

df_pc_train=as.data.frame(predict(pca_fit)) %>% mutate(group=df_train$group)
df_pc_test=as.data.frame(predict(pca_fit,newdata=df_test[features_in])) %>% mutate(group=df_test$group)


#Run KNN ----
print("Run KNN")

formula_in=as.formula(paste0("group ~ ",paste0(paste0("PC",1:pca_dims),collapse = "+")))

knnFit_rest = caret::train(formula_in, data = df_pc_train, method = "knn",
                    preProcess = c("center", "scale"),
                    tuneGrid = expand.grid(k = c(KNN_K)))

df_pc_train$group_pred=predict(knnFit_rest) %>% as.vector()

df_pc_test$group_pred=predict(knnFit_rest,df_pc_test) %>% as.vector()

df_out=data.frame(
  dataset_label=dataset_label,
  seed_in=seed_in,
  n_for_train=n_for_train,
  n_for_features=n_for_features,
  KNN_K=KNN_K,
  pca_dims=pca_dims,
  error_rate_inTrainData=round((nrow(df_pc_train)-sum(df_pc_train$group==df_pc_train$group_pred))/nrow(df_pc_train),8),
  error_rate_inTestData=round((nrow(df_pc_test)-sum(df_pc_test$group==df_pc_test$group_pred))/nrow(df_pc_test),8),
  features_in=paste0(features_in,collapse = ","),
  stringsAsFactors = F
)

write.table(
  data.frame(features_in=features_in) %>% sample_n(10),
  "test/features_in_test.list",
  col.names = F,row.names = F,quote=F
)

#Write output ----
out_file=paste0(dir_out,"GA_",dataset_label,"_seed",seed_in,"_NTrain",n_for_train,"_NFeatures",n_for_features,"_NKNN",KNN_K,'_DimsPca',pca_dims,".tsv")
print(out_file)

write_tsv(df_out,out_file)




