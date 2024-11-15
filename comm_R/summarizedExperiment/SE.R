get_ref_cellType=function(df_ref,key_cell){
  colData_ref=colData(df_ref)
  
  df_label=data.frame(
    label.main=colData_ref$label.main,
    label.fine=colData_ref$label.fine,
    stringsAsFactors = F
  ) %>% mutate(
    cellType=ifelse(label.main==key_cell,label.fine,label.main)
  )
  
  df_label1=df_label["cellType"] %>% distinct()
  df_label1$col=sort(unique(subtypeCol))[1:nrow(df_label1)]
  
  df_label=df_label %>% left_join(df_label1)
  df_label
}












