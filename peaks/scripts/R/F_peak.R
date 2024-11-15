#obj_out=list(SE=SummarizedExperiment(assays=list(vst=df_vst,counts=df_counts), colData = diag))
#assays(obj_in$SE)$counts=count_2042
#assays(obj_in$SE)[['counts']]=count_2042
#obj_in1$SE=obj_in1$SE[,row.names(df_colData)] subset samples
#print(paste0("Input N(samples),N(features): ",paste0(dim(matrix_),collapse = ",")))
#obj_in$SE[["cluster_Phenograph"]]=df_cluster$cluster_Phenograph add new variable

# Import required libraries --------------------------------------------------------------------------
options(warn=-1)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  # library(DESeq2)
  library(stringr)
  library(SummarizedExperiment)
})
options(warn=0)

# install.packages("caret")
# BiocManager::install("SummarizedExperiment")

#args <- commandArgs(trailingOnly = TRUE)
#print(args)

#file_tsv=args[1]

quitely_load_library=function(name){
  suppressMessages(suppressWarnings(library(name,quietly = T,verbose=F,character.only = T)))
}



