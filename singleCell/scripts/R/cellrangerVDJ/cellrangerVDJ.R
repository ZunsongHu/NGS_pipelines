# functions ------------------------------------------------------------------------------------------- ----

get_cdr3_freq=function(object_in,df_vdj_in,top_count_cdr3,label,outdir){
  if(!dir.exists(outdir)){dir.create(outdir)}
  
  
  count_matrix=GetAssayData(object_in)
  count_matrix1=as.data.frame(count_matrix)
  
  df_celltype=object_in@meta.data
  df_celltype$barcode=row.names(df_celltype)
  
  df_celltype1=df_celltype %>% select(barcode,SingleR.labels)
  
  df_cdr3_in1=df_vdj_in %>% left_join(df_celltype1) %>% filter(!is.na(SingleR.labels)) %>%
    filter(productive=="True" & full_length=="True" & is_cell=="True" & high_confidence=="True") %>%
    filter(chain=="IGH") %>% group_by(barcode) %>% arrange(desc(reads)) 
  
  df_freq_cdr3=as.data.frame(table(df_cdr3_in1$cdr3)) %>% arrange(desc(Freq)) %>% 
    transmute(cdr3=as.character(Var1),count_cdr3=Freq)
  
  df_vdj_enriched=df_cdr3_in1 %>% filter(cdr3 %in% head(df_freq_cdr3$cdr3,top_count_cdr3)) %>% 
    left_join(df_freq_cdr3) %>% arrange(desc(count_cdr3),cdr3)
  
  openxlsx::write.xlsx(df_freq_cdr3,file.path(outdir,paste0("df_count_cdr3_",label,".xlsx")),overwrite = T)
  openxlsx::write.xlsx(df_vdj_enriched,file.path(outdir,paste0("df_vdj_enriched_",label,".xlsx")),overwrite = T)
  
  list(df_vdj_enriched=df_vdj_enriched,df_count_cdr3=df_freq_cdr3)
}

#draw enrichment plot ------------------------------------------------------------------------------------------- ----
draw_enrichment_plot=function(indata,label){
  relations1=indata %>% left_join(df_celltype) %>% left_join(df_igh) 
  relations1$SingleR.labels1=paste0(relations1$SingleR.labels,1:nrow(relations1))
  relations1$c_gene1=paste0(relations1$c_gene,1:nrow(relations1))
  relations2=relations1[c("cdr3","c_gene1")]
  
  
  shape1=data.frame(name=unique(relations1$cdr3),shape="circle",size=1,stringsAsFactors = F)
  shape2=relations1[c("c_gene1","shape")] %>% mutate(size=5)
  names(shape2)[1:2]=c("name","shape")
  df_shape=rbind(shape1,shape2)
  
  g <- graph_from_data_frame(relations2,directed=F)
  
  V(g)[relations2$cdr3]$color = "black"
  V(g)[relations2$c_gene1]$color=relations1$col
  
  #pdf(paste0("2.data_temp/vdj_table/",label,"_enrichment_plot.pdf"),width = 8,height = 6)
  png(paste0("2.data_temp/vdj_table/",label,"_enrichment_plot.png"),width = 2400,height = 1800,res = 300)
  plot(g,vertex.label=NA,vertex.size=df_shape$size,vertex.shape=df_shape$shape)
  legend("topleft",pch=19,col=df_celltype$col,legend = df_celltype$celltype1,
         bg=NA,bty = "n")
  legend("left",legend = df_igh1$c_gene1,pch=df_igh1$pch,bg=NA,bty = "n")
  dev.off()
}
