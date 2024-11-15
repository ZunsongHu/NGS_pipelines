#rm(list=ls())
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

# library(SummarizedExperiment, quietly = T)

#testing parameters ----
# setwd("/scratch/zuhu/project/BALL/feature_selection/")
# file_list=system("ls 2042_middleGroup_mRNA/Boruta/*.rds",intern = T)
# label=          "2042_middleGroup_mRNA"
# split_dataset=  "F"
# file_anno=      "/ref_genomes/genome_anno/human/gtf/V102/Homo_sapiens.GRCh38.V102.anno"
# file_featureProperty="2042_middleGroup_mRNA/featureProperty/df_featureProperty.2042_middleGroup_mRNA.tsv"
# dir_out=        "test/"

# "2042_middleGroup_mRNA/Boruta/Boruta_2042_middleGroup_mRNA_EachGroupN100_diag.rds"

#parameters --------------------------------------------------------------------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
print(args)

file_list=            args[1]
label=                args[2]
split_dataset=        args[3]
file_anno=            args[4]
file_featureProperty= args[5]
dir_out=              args[6]

n_digits=4

#create dir out --------------------------------------------------------------------------------------------------------------------------------------------------------
if(!dir.exists(dir_out)){dir.create(dir_out)}

#read boruta output --------------------------------------------------------------------------------------------------------------------------------------------------------
file_list=readLines(file_list)

if(!file_anno=="NULL"){
  df_anno=read_tsv(file_anno)
}

#get featureProperty --------------------------------------------------------------------------------------------------------------------------------------------------------
if(!file_featureProperty=="NULL"){
  df_featureProperty=read_tsv(file_featureProperty)
}


#get boruta summarised table --------------------------------------------------------------------------------------------------------------------------------------------------------

if(!split_dataset=="T"){
  #get all boruta out --------------------------------------------------------------------------------------------------------------------------------------------------------
  df_boruta_out=bind_rows(lapply(file_list, function(x){
    strs=unlist(strsplit(basename(x),"_"))
    EachGroupN=as.numeric(gsub("EachGroupN","",strs[grepl("EachGroupN",strs)]))
    obj_boruta=readRDS(x)
    df_boruta_out=obj_boruta$boruta_out %>% 
      arrange(desc(meanImp)) %>% 
      mutate(
        EachGroupN=EachGroupN,
        rank=1:n())
  }))
  
  #add annotation infor --------------------------------------------------------------------------------------------------------------------------------------------------------
  if("gene_id_raw" %in% names(df_boruta_out)){
    names(df_boruta_out)[match("gene_id_raw",names(df_boruta_out))]="feature_id"} else {
          names(df_boruta_out)[match("gene_id",names(df_boruta_out))]="feature_id"
    }
  
  df_boruta_out$gene_id=NULL
  df_boruta_out$feature_id=gsub("[`]","",df_boruta_out$feature_id)
  df_boruta_out=df_boruta_out %>% left_join(df_anno)

  #convert to digits 
  var_convert=c('meanImp','medianImp','minImp','maxImp','normHits','meanImp_Tentative','medianImp_Tentative','minImp_Tentative','maxImp_Tentative','normHits_Tentative')
  df_boruta_out[var_convert]=apply(df_boruta_out[var_convert], 2, function(x){round(x,n_digits)})
  
  #output raw Boruta output
  cat("Write raw Boruta output\n")
  write_tsv(df_boruta_out,paste0(dir_out,"Boruta_",label,"_rawOut.tsv"))
  
  #Get top genes --------------------------------------------------------------------------------------------------------------------------------------------------------
  # table(df_boruta_out$decision)
  features_pass=unique(df_boruta_out$feature_id[df_boruta_out$decision %in% c('Tentative','Confirmed')])
  
  df_boruta_out1=df_boruta_out %>%
    filter(feature_id %in% features_pass) %>%
    mutate(EachGroupN1=paste0("EachGroupN",EachGroupN)) %>%
    group_by(feature_id) %>%
    mutate(
      rank_mean=round(mean(rank),2),
      rank_median=round(median(rank),2),
      rank_max=round(max(rank),2),
      rank_min=round(min(rank),2),
    ) 
  
  df_boruta_out2=df_boruta_out1[c(names(df_anno),'rank','rank_mean','rank_median','rank_max','rank_min','EachGroupN1')] %>%
    tidyr::spread(EachGroupN1,rank) %>%
    arrange(rank_min,rank_median,rank_mean) 
  
  if(!file_featureProperty=="NULL"){
    df_boruta_out2=df_boruta_out2 %>% left_join(df_featureProperty)
  }

  cat("Write Top list (including Confirmed and tentative ones)\n")
  write_tsv(df_boruta_out2,paste0(dir_out,"Boruta_",label,"_TopList.tsv"))
  
  #Get confirmed ones --------------------------------------------------------------------------------------------------------------------------------------------------------
  features_confirmed=unique(df_boruta_out$feature_id[df_boruta_out$decision %in% c('Confirmed')])
  
  df_confirmed_times=df_boruta_out1[c(names(df_anno),'rank','rank_mean','rank_median','rank_max','rank_min','EachGroupN1','decision')] %>% 
    filter(decision == "Confirmed") %>% 
    group_by(feature_id) %>% 
    mutate(confirmed_times=n()) %>% select(feature_id,confirmed_times) %>% distinct()
  
  df_boruta_confirmed=df_boruta_out1[c(names(df_anno),'rank','rank_mean','rank_median','rank_max','rank_min','EachGroupN1')] %>% 
    filter(feature_id %in% features_confirmed) %>% 
    left_join(df_confirmed_times) %>% 
    tidyr::spread(EachGroupN1,rank) %>%
    arrange(rank_min,rank_median,rank_mean)
  
  if(!file_featureProperty=="NULL"){
    df_boruta_confirmed=df_boruta_confirmed %>% left_join(df_featureProperty)
  }
  
  cat("Write confirmed list\n")
  write_tsv(df_boruta_confirmed,paste0(dir_out,"Boruta_",label,"_ConfirmedList.tsv"))
  
  #Draw confirmed feature heatmap --------------------------------------------------------------------------------------------------------------------------------------------------------
  df_confimed_features=df_boruta_out1 %>% filter(feature_id %in% features_confirmed)
  
  df_confimed_heatmap_plot=df_confimed_features %>% mutate(EachGroupN1=paste0("EachGroupN",sprintf("%03d",EachGroupN))) %>% 
    left_join(df_confirmed_times) %>% 
    select(feature_id,gene_name,EachGroupN1,confirmed_times,rank_median,decision)  %>% 
    spread(EachGroupN1,decision) %>% 
    arrange(rank_median,desc(confirmed_times)) %>% 
    group_by(gene_name) %>% 
    mutate(n_geneName=n(),
           gene_name=ifelse(n_geneName>=2,paste0(gene_name,"_",1:n()),gene_name),
           n_geneName=NULL)  %>% as.data.frame()
  
  write_tsv(df_confimed_heatmap_plot,paste0(dir_out,"Boruta_",label,"_ConfirmedFeaturesHeatmap.tsv"))
  
  row.names(df_confimed_heatmap_plot)=df_confimed_heatmap_plot$gene_name
  
  var_in=c('EachGroupN010','EachGroupN025','EachGroupN050','EachGroupN075','EachGroupN100','EachGroupN150','EachGroupN200','EachGroupN250')
  
  matrix_plot=df_confimed_heatmap_plot[var_in]
  
  gg_dimHeatmap(matrix_plot,y_tick_lab_size = 0.5)
  
  height_=floor(nrow(matrix_plot)*0.01)
  width_=ceiling(height_/2)
  
  ggsave(paste0(dir_out,"Boruta_",label,"_ConfirmedFeaturesHeatmap.pdf"),height = height_,width=width_)
  }


if(split_dataset=="T"){
  df_boruta_out=bind_rows(lapply(file_list, function(x){
    strs=unlist(strsplit(basename(x),"_"))
    split_i=as.numeric(word(x,-2,sep="[/]"))
    EachGroupN=as.numeric(gsub("EachGroupN","",strs[grepl("EachGroupN",strs)]))
    obj_boruta=readRDS(x)
    df_boruta_out=obj_boruta$boruta_out %>% 
      arrange(desc(meanImp)) %>% 
      mutate(
        gene_id=NULL,
        EachGroupN=EachGroupN,
        split_i=split_i,
        rank=1:n())
  }))
  
  write_tsv(df_boruta_out,paste0(dir_out,"Boruta_",label,"_rawOut.tsv"))
  
  df_boruta_out1=df_boruta_out %>% 
    filter(decision %in% c("Confirmed","Tentative")) %>% 
    mutate(gene_id_raw=gsub("[`]","",gene_id_raw)) %>% 
    select(gene_id_raw) %>% 
    distinct()
  
  write.table(df_boruta_out1,paste0(dir_out,"Boruta_",label,"_ConfirmedTentative.list"),col.names = F,row.names = F,quote = F)
  
}
  














