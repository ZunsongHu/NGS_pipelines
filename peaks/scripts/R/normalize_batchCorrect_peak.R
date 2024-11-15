# set ---------------------------
source("/home/zgu_labs/bin/R/functions.R")
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

# testing parameters ----
# setwd("/scratch/zuhu/project/Grp_NadiaCarlesso/atac/")
# 
# dir_out=      "test/"
# batch_var=    "Batch"
# file_obj=     "analysis/second_23samples/peak_counts/obj.second_23samples.filtered.rds"
# parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

dir_out=        args[1]
label=          args[2]
batch_var=      args[3]
file_obj=       args[4]

# create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

# read object ----
obj_=readRDS(file_obj)
df_info=as.data.frame(colData(obj_$SE))

write_tsv(df_info,paste0(dir_out,"df_info.",label,"",".tsv"))


obj_$SE

# run cpm ----
message("Running CPM ...")
logCpm=edgeR::cpm(obj_$SE, lib.size = T, log = T )

# run batch correction ----
if(batch_var=="NULL"){batch_var=NULL}
if(is.null(batch_var)){message("No batch detected. Does not run batch correction!")
  logCpmCombat=logCpm
  }

if(!is.null(batch_var)){if(!batch_var=="NULL"){
  message("Running batch correction ...")
  
  modcombat = model.matrix(~1, data=df_info)
  batch=factor(unlist(df_info[batch_var]))
  
  set.seed(1)
  logCpmCombat = sva::ComBat(dat=as.matrix(logCpm), batch=batch, mod=modcombat, par.prior=TRUE) %>%
    as.data.frame() %>% round(., digits = 4)
}}
  
# get obj ----
message("Get output object ...")
obj_=list(SE=SummarizedExperiment(assays=list(cpm=logCpmCombat), colData = df_info))

obj_$MAD1000=get_variable_genes(obj_in = obj_,N_genes = 1000,assay_name_in = "cpm")

obj_$MAD1000_noCorr=get_nonCorr_genes(obj_in = obj_,feature_panel = "MAD1000",featureRank = 1:1000,highCorrCutoff = 0.75,iteration_times = 2,assay_name_in = "cpm")

saveRDS(obj_,paste0(dir_out,"obj.",label,".peak.cpm",".rds"))

message("\nDone normalization\n")


























