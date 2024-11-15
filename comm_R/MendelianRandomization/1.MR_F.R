########################################### MR function ##################################################

#file1="2.data_temp/out_meta_plink/MotherOverweightObesityVSNormal.out"
#file1_clump="2.data_temp/out_meta_plink/MotherOverweightObesityVSNormal.clumped"
#file2="2.data_temp/out_meta_plink/X1105.out"
#file2_clump="2.data_temp/out_meta_plink/X1105.clumped"
#cutoff=0.00001

#cal_mr_bidirectional(file1,file1_clump,file2,file2_clump,cutoff)
  
read_plink_out_for_MR_filtered=function(file1,file2){
  #Get significant SNPs
  data1=TwoSampleMR::read_exposure_data(filename = file1,snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "A1",pval_col = "P")
  data2=TwoSampleMR::read_exposure_data(filename = file2,snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "A1",pval_col = "P")
  
  
  data1_filtered=TwoSampleMR::read_exposure_data(filename = file1,snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "A1",pval_col = "P")
  data2_filtered=TwoSampleMR::read_exposure_data(filename = file2,snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "A1",pval_col = "P")
  
  #Generate output dataset
  varname1=word(file1,-1,sep="/");varname2=word(file2,-1,sep="/");
  
  data1_1=data1 %>% transmute(snp=SNP,beta_=beta.exposure,se_=se.exposure,p_=pval.exposure)
  data1_filtered1=data1_filtered %>% transmute(snp=SNP,beta_=beta.exposure,se_=se.exposure,p_=pval.exposure)
  
  names(data1_1)[2:4]=paste0(names(data1_1)[2:4],varname1)
  names(data1_filtered1)[2:4]=paste0(names(data1_filtered1)[2:4],varname1)
  
  data2_1=data2 %>% transmute(snp=SNP,beta_=beta.exposure,se_=se.exposure,p_=pval.exposure)
  data2_filtered1=data2_filtered %>% transmute(snp=SNP,beta_=beta.exposure,se_=se.exposure,p_=pval.exposure)
  
  names(data2_1)[2:4]=paste0(names(data2_1)[2:4],varname2)
  names(data2_filtered1)[2:4]=paste0(names(data2_filtered1)[2:4],varname2)
  
  #merge output data
  data_var1_to_var2=data1_filtered1 %>% left_join(data2_1)
  data_var2_to_var1=data2_filtered1 %>% left_join(data1_1)
  
  
  print(paste0('N_snp_var1_to_var2=',nrow(data_var1_to_var2),'; ','N_snp_var2_to_var1=',nrow(data_var2_to_var1)))
  
  return(list(data_var1_to_var2=data_var1_to_var2,data_var2_to_var1=data_var2_to_var1,varname1=varname1,varname2=varname2))
}


read_plink_out_for_MR=function(file1,file1_clump,file2,file2_clump,cutoff){
  #Read full snp data
  data1=TwoSampleMR::read_exposure_data(filename = file1,snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "A1",pval_col = "P")
  data2=TwoSampleMR::read_exposure_data(filename = file2,snp_col = "SNP",beta_col = "BETA",se_col = "SE",effect_allele_col = "A1",pval_col = "P")
  
  #Read independent snp data
  data_clump1=read.table(file1_clump,stringsAsFactors = F,header = T)
  data_clump2=read.table(file2_clump,stringsAsFactors = F,header = T)
  
  #Get significant SNPs
  data1_filtered=data1[data1$pval.exposure<cutoff & data1$SNP %in% data_clump1$SNP, ]
  data2_filtered=data2[data2$pval.exposure<cutoff & data2$SNP %in% data_clump2$SNP & (!data2$SNP %in% data1_filtered$SNP), ]
  
  #Generate output dataset
  varname1=word(file1,-1,sep="/");varname2=word(file2,-1,sep="/");
  
  data1_1=data1 %>% transmute(snp=SNP,beta_=beta.exposure,se_=se.exposure,p_=pval.exposure)
  data1_filtered1=data1_filtered %>% transmute(snp=SNP,beta_=beta.exposure,se_=se.exposure,p_=pval.exposure)
  
  names(data1_1)[2:4]=paste0(names(data1_1)[2:4],varname1)
  names(data1_filtered1)[2:4]=paste0(names(data1_filtered1)[2:4],varname1)
  
  data2_1=data2 %>% transmute(snp=SNP,beta_=beta.exposure,se_=se.exposure,p_=pval.exposure)
  data2_filtered1=data2_filtered %>% transmute(snp=SNP,beta_=beta.exposure,se_=se.exposure,p_=pval.exposure)
  
  names(data2_1)[2:4]=paste0(names(data2_1)[2:4],varname2)
  names(data2_filtered1)[2:4]=paste0(names(data2_filtered1)[2:4],varname2)
  
  #merge output data
  data_var1_to_var2=data1_filtered1 %>% left_join(data2_1)
  data_var2_to_var1=data2_filtered1 %>% left_join(data1_1)
  
  
  print(paste0('N_snp_var1_to_var2=',nrow(data_var1_to_var2),'; ','N_snp_var2_to_var1=',nrow(data_var2_to_var1)))
  
  return(list(data_var1_to_var2=data_var1_to_var2,data_var2_to_var1=data_var2_to_var1,varname1=varname1,varname2=varname2))
}


cal_mr_oneDirection=function(indata,var1,var2){
  n_snp=nrow(indata)
  
  beta_x=unlist(indata[paste0("beta_",var1)])
  se_x=unlist(indata[paste0("se_",var1)])
  
  beta_y=unlist(indata[paste0("beta_",var2)])
  se_y=unlist(indata[paste0("se_",var2)])
  
  sig_snp=paste0(indata$snp,collapse = ";")
  
  if(n_snp>=3){
    MRInputObject <- MendelianRandomization::mr_input(bx = beta_x,
                                                      bxse= se_x,
                                                      
                                                      by = beta_y,
                                                      byse = se_y)
    
    
    MRAllObject_main <- MendelianRandomization::mr_allmethods(MRInputObject, method = "main")
    
    
    mr_out=MRAllObject_main@Values %>% 
      mutate(Beta_se=paste0(round(Estimate,6)," (",round(`Std Error`,6),")"),P=signif(`P-value`,2),
             method=gsub("-","",gsub(" ","",Method)),newvar="x")
    mr_out=mr_out %>% mutate(method=gsub("[()]","",method)) %>% 
      filter(method %in% c("Simplemedian" , "Weightedmedian","IVW","MREgger", "intercept"))
    mr_out1=dcast(mr_out,newvar~method,value.var = "Beta_se")
    mr_out2=dcast(mr_out,newvar~method,value.var = "P")
    mr_out1$newvar=NULL;mr_out2$newvar=NULL
    names(mr_out1)=paste0(names(mr_out1),"_BetaSe");names(mr_out2)=paste0(names(mr_out2),"_P")
    
    
    out=cbind(var1,var2,n_snp,mr_out1,mr_out2)
    out$var1=as.character(out$var1);out$var2=as.character(out$var2);out$n_snp=as.character(out$n_snp)
  }
  
  if(n_snp<=2){
    out=data.frame(var1=var1,var2=var2,n_snp=n_snp,stringsAsFactors = F)
  }
  out$sig_snp=sig_snp
  out$n_snp=as.numeric(out$n_snp)
  return(out)
}
 


#plink_out_pair_one=read_plink_out_for_MR(file1,file1_clump,file2,file2_clump,cutoff)

#indata=plink_out_pair_one$data_var1_to_var2
#indata=plink_out_pair_one$data_var2_to_var1

#var1=plink_out_pair_one$varname1
#var2=plink_out_pair_one$varname2

  

cal_mr_biDirection=function(file1,file1_clump,file2,file2_clump,cutoff){
  plink_out_pair_one=read_plink_out_for_MR(file1,file1_clump,file2,file2_clump,cutoff)
  bind_rows(
    cal_mr_oneDirection(plink_out_pair_one$data_var1_to_var2,plink_out_pair_one$varname1,plink_out_pair_one$varname2),
    cal_mr_oneDirection(plink_out_pair_one$data_var2_to_var1,plink_out_pair_one$varname2,plink_out_pair_one$varname1)
  )
}

#cal_mr_biDirection(file1,file1_clump,file2,file2_clump,cutoff)


#save(read_plink_out_for_MR,cal_mr_oneDirection,cal_mr_biDirection,file = "2.data_temp/mr_function.rdata")



########################################### Code on server  ##################################################
#load("mr_function.rdata")

#library(reshape2)
#library(stringr)
#library(dplyr)

#args = commandArgs(trailingOnly=TRUE)

#var1=args[1]
#var2=args[2]
#cutoff=as.numeric(args[3])

#file1=paste0("3.out_plink_raw/",var1,".out")
#file1_clump=paste0("4.out_clump/",var1,".clumped")

#file2=paste0("3.out_plink_raw/",var2,".out")
#file2_clump=paste0("4.out_clump/",var2,".clumped")

#out=cal_mr_biDirection(file1,file1_clump,file2,file2_clump,cutoff)
#write.csv(out,paste0('5.out_mr_bidirectional/',var1,'_and_',var2,"_cutoff",cutoff,'.csv'),quote=F,row.names=F)













































