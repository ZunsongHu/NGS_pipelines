#get_sampleRelatedness_label ---------------

get_sampleRelatedness_label=function(indf){
  indf1=indf %>% transmute(COH_sample=ID1,InfType=paste0(ID2,"(",InfType,")"))
  indf2=indf %>% transmute(COH_sample=ID2,InfType=paste0(ID1,"(",InfType,")"))
  
  out_df=bind_rows(indf1,indf2) %>% dplyr::mutate(varname="King_relate") %>% 
    reshape2::dcast(COH_sample~varname,function(x){paste0(x,collapse = ",")},value.var = "InfType") 
  out_df
}
