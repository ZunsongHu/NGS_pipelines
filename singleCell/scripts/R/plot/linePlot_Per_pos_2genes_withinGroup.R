
gene="Cd69"

w_f=5
h_f=4.5



#Cd69 -------------------------------
x="per_prlr"
y="per_gene"

obj_logC$celltype_mergeB=ifelse(as.character(obj_logC$celltype_manual) %in% c('Naive B cell','Memory B cell',"Plasma cell"),"B cell",as.character(obj_logC$celltype_manual))

labels_group=c("Health" = "Healthy spleen (n=3)", "Control" = "SLE spleen Control SMO (n=4)", "Antimia" = "SLE spleen LFPRLR SMO (n=4)")


get_linPlot=function(obj_in=obj_logC,gene="Cd69",g_celltype="celltype_mergeB",
                     x_lab="Prop of Prlr+ cells",y_lab="Prop of Cd69+ cells",
                     
                     dir_out="analysis/linePlot/Prlr_Cd69/"){
  
  cat(paste0("Working for ",gene))
  
  df_count=FetchData(obj_in,assay="RNA",layer="counts",vars = c("group","group1","group2","Prlr_exp",g_celltype,gene))
  names(df_count)[c(ncol(df_count)-1,ncol(df_count))]=c("g_celltype","gene")
  
  
  df_celltype1= df_count %>% 
    mutate(
      gene_exp=ifelse(gene>0,"Exp","noExp"),
    ) %>% 
    filter(!g_celltype=="Erythrocyte") %>% 
    group_by(group,g_celltype) %>% mutate(n=n())
  
  
  df_g=df_celltype1 %>% dplyr::select(group,group1,g_celltype) %>% distinct()
  
  df_prlr=df_celltype1 %>% group_by(group,g_celltype,Prlr_exp) %>% mutate(n_g=n(),per_prlr=n_g/n) %>% 
    dplyr::select(group,group1,g_celltype,Prlr_exp,per_prlr) %>% distinct() %>% filter(Prlr_exp=="Exp")
  
  df_gene=df_celltype1 %>% group_by(group,g_celltype,gene_exp) %>% mutate(n_g=n(),per_gene=n_g/n) %>% 
    dplyr::select(group,group1,g_celltype,gene_exp,per_gene) %>% distinct() %>% filter(gene_exp=="Exp")
  
  df_gene_per=df_g %>% left_join(df_prlr) %>% left_join(df_gene) %>% 
    mutate(
      per_prlr=ifelse(is.na(Prlr_exp),0,per_prlr),
      per_gene=ifelse(is.na(per_gene),0,per_gene),
    )
  
  celltypes=sort(unique(df_gene_per$g_celltype))
  celltype=celltypes[1]
  celltype="Plasma cell"
  
  for(celltype in celltypes){
    
    df_in=df_gene_per[df_gene_per$g_celltype == celltype,]
    df_in=df_in[c("per_prlr","per_gene","group1")]
    names(df_in)=c('x','y',"group1")
    
    fit_lm = lm(y ~ x, df_in)
    fit_results=summary(fit_lm)
    df_p=fit_results$coefficients %>% as.data.frame()
    
    range_x=(max(df_in$x)-min(df_in$x))/5
    
    P=paste0("P = ",signif(df_p$`Pr(>|t|)`[2],2))
    
    title_=paste0(celltype," w.r.t total WBC\n",P)
    
    ggplot(df_in,aes(x=x,y=y,color=group1))+
      geom_point()+
      labs(title = title_,x=x_lab,y=y_lab,color=NULL)+
      scale_color_manual(values=col_g1,labels=labels_group)+
      theme_paper+
      theme(plot.title = element_text(hjust = 0.5, size = 15),legend.position = "bottom")+
      guides(color = guide_legend(ncol = 1))+
      geom_abline(slope = coef(fit_lm)[2], 
                  intercept = coef(fit_lm)[["(Intercept)"]])+
      legend_axis_text_size(12)

    ggsave(paste0(dir_out,celltype,".png"),width = w_f,height=h_f)
    ggsave(paste0(dir_out,celltype,".pdf"),width = w_f,height=h_f)
  }
}

get_linPlot(obj_in=obj_logC,gene="Cd69",g_celltype="celltype_mergeB",
                     x_lab="Prop of Prlr+ cells w.r.t\n total WBC",y_lab="Prop of Cd69+ cells w.r.t\n total WBC",
                     dir_out="analysis/linePlot/Prlr_Cd69_celltypeB/")

get_linPlot(obj_in=obj_logC,gene="Cd69",g_celltype="celltype_manual",
                     x_lab="Prop of Prlr+ cells w.r.t\n total WBC",y_lab="Prop of Cd69+ cells w.r.t\n total WBC",
                     dir_out="analysis/linePlot/Prlr_Cd69/")