options(future.globals.maxSize =10*1024 * 1024^2) #10GB
gc(rm(list=ls()))
rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

suppressWarnings(suppressMessages(require(SummarizedExperiment)))
suppressWarnings(suppressMessages(require(MLSeq)))
suppressWarnings(suppressMessages(require(DaMiRseq)))
suppressWarnings(suppressMessages(require(GARS)))

#parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

obj_rds=          args[1]
dataset_label=    args[2]
var_group=        args[3]
dir_out=          args[4]

#testing parameters ----
setwd("/scratch/zuhu/project/BALL/feature_selection/")

obj_rds=          "2042/smote/smote_2042_diag_EachGroupN10.rds"
dataset_label =   "2042_EachGroupN50"
var_group=        "diag"
dir_out=          "test/"


obj_type=         "matrix"
features=gene_in
# features=features[1:100]
#create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read obj ----
obj_=readRDS(obj_rds)

if(tolower(obj_type)=="matrix"){
  value_group=unlist(obj_[var_group])
  df_in=obj_[c(features)]
  # df_in$group=value_group
  }

#Run GARS ----
print("Run GARS")
set.seed(0)

res_GA <- GARS_GA(data=df_in,
                  classes = value_group,
                  chr.num = 1000,
                  chr.len = 100,
                  generat = 20,
                  co.rate = 0.8,
                  mut.rate = 0.1,
                  n.elit = 10,
                  type.sel = "RW",
                  type.co = "one.p",
                  type.one.p.co = "II.quart",
                  n.gen.conv = 150,
                  plots="no",
                  verbose="yes")


data_reduced_GARS = MatrixFeatures(res_GA)

#get output -----
boruta_out_Tentative = TentativeRoughFix(boruta_train)

importance_boruta=attStats(boruta_train)
importance_boruta$gene_id=row.names(importance_boruta)

importance_boruta_Tentative=attStats(boruta_out_Tentative)
names(importance_boruta_Tentative)=paste0(names(importance_boruta_Tentative),"_Tentative")
importance_boruta_Tentative$gene_id=row.names(importance_boruta_Tentative)


df_out=info_gtf %>% 
  left_join(importance_boruta) %>% 
  left_join(importance_boruta_Tentative) %>%
  filter(!is.na(meanImp))
  

#save output ----
obj_out=list(fit_boruta=boruta_train,
             boruta_out=df_out) 


print("Save to output")
saveRDS(obj_out,paste0(dir_out,"Boruta_",dataset_label,"_",var_group,".rds"))



#SMOTE oversampling of minor subtypes
# smoteToNums=c(10, 25, 50, 75, 100, 125, 150, 200)

# library(UBL)
# 
# install.packages("UBL")
# 
# 
# library(DMwR)
# 
# smoteToNums=c(200)
# for (smoteToNum in smoteToNums) {
# }
# 
# smoteID=paste0("Smote-", smoteToNum)
# borutaID=paste0("Smote-", smoteToNum, "_Boruta")
# print(smoteID)
# 







