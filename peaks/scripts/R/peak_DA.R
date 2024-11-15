# set ---------------------------
source("/home/zgu_labs/bin/R/functions.R")
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

library(edgeR)
# testing parameters ----
setwd("/scratch/zuhu/project/Grp_NadiaCarlesso/atac/")

dir_out=      "test/"
label=        "test"
batch_var=    "Batch"
compare_var=  "group"
file_obj=     "analysis/second_23samples/peak_counts/obj.second_23samples.filtered.rds"
file_counts=  "analysis/peak_counts/peakCounts.All.all.rds"
# parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

dir_out=        args[1]
label=          args[2]
batch_var=      args[3]
compare_var=    args[4]
file_obj=       args[5]

# create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

# read object ----
obj_filtered=readRDS(file_obj)
df_info=as.data.frame(colData(obj_filtered$SE))

# read counts ----
counts=readRDS(file_counts)

#get analyzing object ----
df_row=rowRanges(counts) %>% as.data.frame() %>% mutate(feature=paste0(seqnames,":",start,"-",end))
df_col=counts@colData %>% as.data.frame() %>% mutate(file=basename(bam.files),COH_sample=word(file,1,sep="[.]")) 

df_counts=assays(counts)[["counts"]]
row.names(df_counts)=df_row$feature
colnames(df_counts)=df_col$COH_sample

df_counts_=df_counts[,df_info$COH_sample]
dim(df_counts_)
obj_=list(SE=SummarizedExperiment(assays=list(counts=df_counts_), colData = df_info))


# get compare variable level ----
df_info$compare_var=unlist(df_info[compare_var])
if(any(is.na(df_info$compare_var))){stop("Missing values found in group for comparison")}

g_n=df_info %>% group_by(compare_var) %>% mutate(n=n()) %>% select(compare_var,n) %>% distinct()
if(any(g_n$n<2)){stop("All freqs of each must more than or equal to 2")}

compare_var_levels=unique(g_n$compare_var)

# one level vs all the rest as whole ----
message("\nRunning one level vs all the rest as whole ...\n")

dir.create(paste0(dir_out,"oneLevel_vs_others/"))

# compare_var_level=compare_var_levels[1]

for(compare_var_level in compare_var_levels){
  df_DE_=run_DE_edgeR(obj_in = obj_,var_effect = compare_var,test_level = compare_var_level,control_level = "NULL",var_adj = "Batch") %>% 
    mutate(logFdrP=-log(P_fdr,base = 10),logP=-log(PValue,base = 10))
  write_tsv(df_DE_,paste0(dir_out,"oneLevel_vs_others/df_DA.",label,".",compare_var_level,"_vs_others",".tsv"))
  
  gg_geomPoint(df_in = df_DE_,x = "logFC",y = "logFdrP",size = 1,
                  x_refLine = ,y_refLine = -log(0.05,10),x_lab = "logFC",
                  y_lab = expression(paste("-", log[10]," Adjusted P-value")),
                  title=paste0(compare_var_level,"_vs_others, ","P(fdr)<0.05 N= ",nrow(df_DE_[df_DE_$P_fdr<0.05,])))
  ggsave(paste0(dir_out,"oneLevel_vs_others/volcanoPlot.Pfdr.",label,".",compare_var_level,"_vs_others",".png"),width=6,height = 5)
  ggsave(paste0(dir_out,"oneLevel_vs_others/volcanoPlot.Pfdr.",label,".",compare_var_level,"_vs_others",".pdf"),width=6,height = 5)
  
  gg_geomPoint(df_in = df_DE_,x = "logFC",y = "logP",size = 1,
               x_refLine = ,y_refLine = -log(0.05,10),x_lab = "logFC",
               y_lab = expression(paste("-", log[10]," P-value")),
               title=paste0(compare_var_level,"_vs_others, ","P <0.05 N= ",nrow(df_DE_[df_DE_$PValue<0.05,])))
  ggsave(paste0(dir_out,"oneLevel_vs_others/volcanoPlot.P.",label,".",compare_var_level,"_vs_others",".png"),width=6,height = 5)
  ggsave(paste0(dir_out,"oneLevel_vs_others/volcanoPlot.P.",label,".",compare_var_level,"_vs_others",".pdf"),width=6,height = 5)
  
}

# one level vs all the rest individually ----
message("\nRunning one level vs rest individually ...\n")

dir.create(paste0(dir_out,"pairwise/"))

for(compare_var_level1 in compare_var_levels){
  for(compare_var_level2 in compare_var_levels){
    if (compare_var_level1==compare_var_level2){break}
    df_DE_=run_DE_edgeR(obj_in = obj_,var_effect = compare_var,test_level = compare_var_level1,control_level = compare_var_level2,var_adj = "Batch") %>% 
      mutate(logFdrP=-log(P_fdr,base = 10),logP=-log(PValue,base = 10))
    write_tsv(df_DE_,paste0(dir_out,"pairwise/df_DA.",label,".",compare_var_level1,"_vs_",compare_var_level2,".tsv"))
    
    gg_geomPoint(df_in = df_DE_,x = "logFC",y = "logFdrP",size = 1,
                 x_refLine = ,y_refLine = -log(0.05,10),x_lab = "logFC",
                 y_lab = expression(paste("-", log[10]," Adjusted P-value")),
                 title=paste0(compare_var_level1,"_vs_",compare_var_level2,", P(fdr)<0.05 N= ",nrow(df_DE_[df_DE_$P_fdr<0.05,])))
    ggsave(paste0(dir_out,"pairwise/volcanoPlot.Pfdr.",label,".",compare_var_level1,"_vs_",compare_var_level2,".png"),width=6,height = 5)
    ggsave(paste0(dir_out,"pairwise/volcanoPlot.Pfdr",label,".",compare_var_level1,"_vs_",compare_var_level2,".pdf"),width=6,height = 5)
    
    gg_geomPoint(df_in = df_DE_,x = "logFC",y = "logP",size = 1,
                 x_refLine = ,y_refLine = -log(0.05,10),x_lab = "logFC",
                 y_lab = expression(paste("-", log[10]," P-value")),
                 title=paste0(compare_var_level1,"_vs_",compare_var_level2,", P <0.05 N= ",nrow(df_DE_[df_DE_$PValue<0.05,])))
    ggsave(paste0(dir_out,"pairwise/volcanoPlot.P.",label,".",compare_var_level1,"_vs_",compare_var_level2,".png"),width=6,height = 5)
    ggsave(paste0(dir_out,"pairwise/volcanoPlot.P.",label,".",compare_var_level1,"_vs_",compare_var_level2,".pdf"),width=6,height = 5)
  }}

message("\nDone DA\n")



















