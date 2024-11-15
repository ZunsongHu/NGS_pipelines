




# functions ------------------------------------------------------------------------------------------- ----
draw_DimPlot=function(inObj,df_feature,var_key,label,file_out_prefix){
  if (!all(row.names(inObj@meta.data)==unlist(df_feature[1]))){stop("barcode not equal")}
  
  inObj[[label]]=unlist(df_feature[var_key])
  DimPlot(inObj,group.by = label,label=T,repel = T,cols=DiscretePalette(length(unique(unlist(df_feature[var_key])))))
  
  ggsave(paste0(file_out_prefix,".png"),width=10,height = 10)
  ggsave(paste0(file_out_prefix,".pdf"),width=10,height = 10)
  inObj
}

draw_umap_mini=function(obj,var_group,label,dir_out,width=15,height=11){
  
  n_levels=length(unique(unlist(get_metadataVar(obj,var_group)[,2])))
  
  DimPlot(obj,group.by = var_group,label=T,repel = T,cols=DiscretePalette(n_levels))
  ggsave(paste0(dir_out,"umap_",label,"_",var_group,".pdf"),width=width,height = height)
  ggsave(paste0(dir_out,"umap_",label,"_",var_group,".png"),width=width,height = height)
  
}


draw_umap_class=function(inObj,var,outfile_prefix="2.data_temp/test/VlnPlot_qc"){
  dir_name=dirname(dirname(outfile_prefix));if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  dir_name=dirname(outfile_prefix);if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  var_value=unlist(inObj@meta.data[var])
  
  Idents(inObj)=var_value
  
  col_in=unique(subtypeCol)[2:length(unique(subtypeCol))][1:length(unique(var_value))]
  
  DimPlot(inObj,cols = col_in)
  ggsave(filename = paste0(outfile_prefix,".png"),width = 13,height = 9)
  ggsave(filename = paste0(outfile_prefix,".pdf"),width = 13,height = 9)
}



draw_feature_plot_fromDF=function(objIn,df,varname,label,cols_in,dir_out){
  
  objIn[[label]]=unlist(df[varname])
  
  DimPlot(objIn,group.by = label,label=T,repel = T,
          cols = cols_in[names(cols_in) %in% unlist(objIn@meta.data[label])])
  ggsave(paste0(dir_out,label,".png"),width=16,height = 9)
  ggsave(paste0(dir_out,label,".pdf"),width=16,height = 9)
  
  objIn
}



