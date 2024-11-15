
#get_design_MeRIPseqPipe ---------------------------------------------------------------------------
get_design_MeRIPseqPipe=function(file_merge,file_out){
  
  df_merge=read_tsv(file_merge) %>% 
    dplyr::select(file,read,Sample_ID) %>% 
    mutate(
      merip_g=word(Sample_ID,2,sep="[_-]"),
      merip_g=paste0(ifelse(merip_g=="Input","input",ifelse(merip_g=="caR","ip",NA)),"_R",read),
      Group_ID=word(Sample_ID,2,sep="-"),
      Sample_ID=paste0(Group_ID,word(Sample_ID,-1,sep="-"))
    ) %>% 
    tidyr::pivot_wider(id_cols=c("Sample_ID","Group_ID"), names_from = "merip_g",values_from = "file") %>% 
    select(Sample_ID,input_R1,input_R2,ip_R1,ip_R2,Group_ID)
  
  write_tsv(df_merge,file_out)
}
