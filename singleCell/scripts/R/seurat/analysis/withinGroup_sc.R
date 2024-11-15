#usage:
# R_seurat5 /net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/singleCell/scripts/R/seurat/analysis/withinGroup_sc.R \

options(warn=-1)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(stringr)
  library(scCustomize)

})
options(warn=0)


# functions ------------------------------------------------------------------------------------------------ ----


#parameters ------------------------------------------------------------------------------------------------ ----
setwd("/net/nfs-irwrsrchnas01/labs/zgu_grp/Group/Grp_Swaminathan/singleCell")
file_obj="/net/nfs-irwrsrchnas01/labs/zgu_grp/Group/Grp_Swaminathan/singleCell/analysis/obj/obj_subset/obj_prlrPos.antimiaControl.1491cells.rds"

var_group="celltype_manual"
cols_group=cols

var_compare="group1"
cols_comapre=col_g1

features=c("Stat1","Stat2","Irf7","Oas1a","Oas2","Oas3","Oasl1","Ifitm3","Ifi35","Ifi44","Ifi44l","Ifit1","Ifit3")

genes_boxPlot=c("Stat1","Stat2","Irf7","Oas1a","Oas2","Oas3","Oasl1","Ifitm3","Ifi35","Ifi44","Ifi44l","Ifit1","Ifit3")
prefix=paste0("analysis/withinGrorup/prlrPos.",var_group,".",var_compare,"/")



dir.create(prefix)

# DE
run_DE=F

# boxPlot
run_boxPlot=T


# vlnPlot
run_vln=T
w_vln=20;h_vln=15

obj_$celltype_manual


subset_n=200




#read obj ------------------------------------------------------------------------------------------------ ----

obj_=readRDS(file_obj)

df_bc=get_barcode_df(obj_)

if(!is.null(subset_n)){
  set.seed(10)
  df_bc1=df_bc %>% sample_n(subset_n)
  obj_=subset(obj_,cells=df_bc1$barcode)
}

# get df  ------------------------------------------------------------------------------------------------ ----

df_=FetchData(object = obj_,vars=c(var_group,var_compare,genes_boxPlot))
names(df_)[1:2]=c("var_group","var_compare")

obj_$var_group=df_$var_group
obj_$var_compare=df_$var_compare

table(df_$var_group,df_$var_compare)
levels_=as.character(unique(df_$var_group))


level_1="T cell"
level_1=levels_[1]

# vlnPlot  ------------------------------------------------------------------------------------------------ ----
if(run_vln){
  dir_vln=paste0(prefix,"vlnPLot/")
  dir.create(dir_vln)
  VlnPlot(obj_,group.by = "var_group",split.by =  "var_compare",features = genes_boxPlot,cols = cols_comapre)
  ggsave(paste0(dir_vln,"vlnPlot.gene.png"),w=w_vln,h=h_vln)
  ggsave(paste0(dir_vln,"vlnPlot.gene.pdf"),w=w_vln,h=h_vln)
}



run_analysis_oneLevel=function(level_1){
  print(level_1)
  
  obj_subset=subset(obj_,subset=var_group==level_1)
  obj_subset$var_compare=as.character(obj_subset$var_compare)
  df_freq_x=table(obj_subset$var_compare) %>% as.data.frame()
  
  # print(df_freq_x)
  
  #DE ----
  if(run_DE){
    if(min(df_freq_x$Freq)>=3 & nrow(df_freq_x)>=2){
      Idents(obj_subset)="var_compare"
      # obj_subset=PrepSCTFindMarkers(obj_subset)
      df_DE_=FindMarkers(obj_subset,ident.1 = "Antimia",ident.2 = "Control",min.pct=0,min.diff.pct=0)
      df_DE_1=df_DE_ %>% tibble::rownames_to_column(var = 'feature') %>%  mutate(contrast="Antimia_vs_Control")
      write_tsv(df_DE_1,paste0(prefix,"df_DE.in_",level_1,".tsv"))
      
      gg_volcanoPlot(df_in = df_DE_1,P = "p_val",FC = "avg_log2FC",var_label = "feature",P_cutoff = 1,n_top_to_show = 20,title = paste0("In ",level_1),y_lab = "nominal P value",
                     prefix = paste0(prefix,"df_DE.in_",level_1),w = 8,h = 6)
    }
  }

  
  #boxPlot ----
  if(run_boxPlot){
    df_1=df_ %>% filter(var_group==level_1)
    
    for (feature in features){
      
      gg_boxPlot_dodge_withP(df_in = df_1,var_value=feature,var_group="var_compare",var_group_color=cols_comapre,cols_color=NULL,x_lab=NULL,y_lab=NULL,plot_title=NULL,
                                      x_tick_label_var=NULL,fontSize=14,x_rotate=65,legend_pos="bottom",legend_label=NULL,legend_title=NULL,
                                      add_P=T,add_points=T,
                                      prefix=NULL,w=10,h=6,sig_bar_move=0.05)
    }
    

  }

}

df_1$var_compare


feature=features[1]


















































