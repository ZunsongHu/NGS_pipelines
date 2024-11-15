


source("/net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/comm_R/3.plot/ggplot2/gg.R")
source("/net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/comm_R/3.plot/color_R.R")


get_barPlot_geneExp=function(obj_in,group,gene,cutoff=0,prefix=NULL){
  df_=FetchData(obj_in,vars = c(group,gene),layer = "counts") %>% 
    tibble::rownames_to_column(var = "barcode")
  names(df_)[2:ncol(df_)]=c("group","gene")
  
  df_1=df_ %>% mutate(g_exp=ifelse(gene>cutoff,"Exp","NoExp"))
  df_per=cal_freqPer_2groups(df_in = df_1,level_low = "g_exp",level_high = "group")
  g1=gg_barplot_stack(df_in = df_per,var_group_x = "group",var_group_bar = "g_exp",var_barLabel = "per2",x_rotate = 45)
  
  if(!is.null(prefix)){
    ggsave(paste0(prefix,".boxPlot.png"),plot =g1,width = w,height = h)
    ggsave(paste0(prefix,".boxPlot.pdf"),plot =g1,width = w,height = h)
  }
  
  g1
}

get_barplot_stack=function(obj_in,var_group_x,var_group_bar,prefix=NULL,w=5.5,h=5){
  df_in=get_metadataVar(obj_in,c(var_group_x,var_group_bar)) 
  
  df_freq=cal_freqPer_2groups(df_in,level_small = var_group_bar,level_big = var_group_x)
  
  g1=gg_barplot_stack(df_in = df_freq,var_group_x = var_group_x,var_group_bar = var_group_bar,var_barLabel = "per2",
                   cols_fill = cols_B,x_tick_label_var = labels_group,x_rotate = 45)
  
  if(!is.null(prefix)){
    ggsave(paste0(prefix,".png"),width = w,height = h)
    ggsave(paste0(prefix,".pdf"),width = w,height = h)
    write_tsv(df_freq,paste0(prefix,".tsv"))
  }
  return(list(df_freq,g1))
}





























