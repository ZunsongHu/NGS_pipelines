#get bed file from feature ------------------------------------------------------------------------- ----
get_bedFile=function(indata,var_name="feature",file_out){
  # names(indata)[1]="feature"
  
  indata1=indata[,var_name]
  
  indata1=indata %>% 
    mutate(
      chr=word(feature,1,sep=":"),
      x=word(feature,2,sep=":"),
      start=as.numeric(word(x,1,sep="-"))-1,
      end=word(x,2,sep="-")
    )
  
  write.table(indata1[c('chr','start','end')],file_out,col.names = F,row.names = F,quote = F,sep = "\t")
}
