#Imputation with multi-imputation
data_cell[1:23]=mice::complete(mice::mice(data_cell[1:23],m=1))


######################################## Imputation ##################################
imputation_knn=function(indata){
  set.seed(1)
  row.names(indata)=indata$studyid
  indata["studyid"]=NULL
  indata=t(indata)
  outdata=impute::impute.knn(indata)$data
  outdata=t(outdata)
  studyid=row.names(outdata)
  outdata1=as.data.frame(cbind(studyid,outdata))
  write.csv(outdata1,paste0(outdir,"data_for_filter.csv"),row.names=F)
  outdata1=read.csv(paste0(outdir,"data_for_filter.csv"))
  outdata1
}

imputation_RandomForest=function(indata){
  row.names(indata)=indata$studyid
  indata["studyid"]=NULL
  set.seed(1)
  cl <- makeCluster(3)  
  registerDoParallel(cl) 
  outdata=missForest::missForest(xmis = indata,verbose = TRUE,parallelize = "forests")$ximp
  stopCluster(cl)
  
  studyid=row.names(outdata)
  outdata1=as.data.frame(cbind(studyid,outdata))
  outdata1
}


######################################## Data Filtering ##################################
cal_rsd=function(x){
  round(sd(x,na.rm=T)/mean(x,na.rm=T),3)}

filter_rsd=function(indata,rsd_cutoff,rsd_label){
  rsd_value=data.frame(
    names(indata)[-1],
    apply(indata[,-1],2,cal_rsd),stringsAsFactors = F
  )
  names(rsd_value)=c("Comp_id1",paste0("rsd_after_imputation_",rsd_label))
  data_annotation=data_raw_annotation %>% left_join(rsd_value)
  xlsx::write.xlsx(as.data.frame(data_annotation),paste0(outdir,"data_annotation_",rsd_label,".xls"),row.names = F,sheetName = "1",showNA = F)
  var_keep=rsd_value$Comp_id1[!rsd_value$rsd_after_imputation<rsd_cutoff]
  outdata=indata[c("studyid",var_keep)]
  outdata
}