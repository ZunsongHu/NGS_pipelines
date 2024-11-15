########################################  Convert NA to Zero ########################################
convert_NA_to_zero=function(x){
  ifelse(is.na(x),0,x)
}




########################################  Convert zero to NAs ########################################

convert_zero_to_NA=function(x){
  ifelse(x==0,NA,x)
}

########################################  Convert values to labels ########################################

convert_value_to_label=function(x,value,label){
  for(i in value){
    x1=rep(NA,length(x))
    for(i in value){
      x1[x==i]=label[match(i,value)]
    }
    x1
  }
  x1
}


########################################  paste with removing NAs and sort the values ########################################

paste_remove_NA_with_sort=function(x,link){
  
  if(any(!is.na(x))|
     any(x=="",na.rm = T)
  ){
    x1=x[!is.na(x)]
    x1=x[!x==""]
    x1=sort(x1)
    x2=paste0(x1,collapse = link)
  }
  
  if(all(is.na(x))|all(x=="")){x2=NA}
  x2
}

######################################## Recode value based on zero and median value:0 1 2 ##################################
#x=data_pheno_withCellType$freq_lcaGroup2
convert_ZeroMedian_012=function(x){
  ifelse(x==0,0,
         ifelse(x>median(x[x>0],na.rm = T),2,
                ifelse(x<=median(x[x>0],na.rm = T),1,NA)))
}


######################################## Recode data based on zero and median value:0 1 2 ##################################
#NA will be trasformed as NA
data_convert_ZeroMedian_012=function(indata){
  out=indata
  out[2:ncol(out)]=
    apply(out[2:ncol(out)],2,function(x){
      ifelse(x==0,0,
             ifelse(x>=median(x[!x==0],na.rm=T),2,
                    ifelse(x<median(x[!x==0],na.rm=T),1,NA)))
    })
  out
}

######################################## Recode data based on Missing status and median value:0 1 2 ##################################

data_convert_missingStatus_012=function(indata){
  out=indata
  out[2:ncol(out)]=
    apply(out[2:ncol(out)],2,function(x){
      ifelse(is.na(x),0,
             ifelse(x>=median(x,na.rm=T),2,
                    ifelse(x<median(x,na.rm=T),1,NA)))
    })
  out
}

vector_convert_missingStatus_012=function(x){
  ifelse(is.na(x),0,
         ifelse(x>=median(x,na.rm=T),2,
                ifelse(x<median(x,na.rm=T),1,NA)))
}


######################################## cal dataset missing rate for variables in dataset ##################################
cal_missing=function(indata,var_name_one){
  N_total=nrow(indata)
  N_missing=N_total-table(is.na(unlist(indata[var_name_one])))[1]
  data.frame(
    var_name=var_name_one,
    N_total=nrow(indata),
    N_missing=N_total-table(is.na(unlist(indata[var_name_one])))[1],
    Missing_rate=paste0(sprintf(paste0("%.",2,"f"),N_missing/N_total*100),"%"),
    stringsAsFactors = F
  )
}

######################################## regroup_median mean ##################################

regroup_median_value=function(x){
  ifelse(x>=median(x,na.rm = T),paste0(1," >=",median(x,na.rm = T)),paste0(0," <",median(x,na.rm = T)))
}


regroup_median_level=function(x){
  ifelse(x>median(x,na.rm = T),"High","Low")
}



regroup_mean_level=function(x){
  ifelse(x>=mean(x,na.rm = T),"High","Low")
}



regroup_values=function(x,cut_c){
  min_value=min(x,na.rm=T)
  max_value=max(x,na.rm = T)
  
  if(min_value< min(cut_c)){cut_c1=c(min_value,cut_c)}
  if(min_value>= min(cut_c)){cut_c1=c(min_value,cut_c[-1])}
  
  if(max_value< max(cut_c)){cut_c2=c(cut_c1[-length(cut_c1)],max_value)}
  if(max_value>= max(cut_c)){cut_c2=c(cut_c1,max_value)}
  
  
  x1=as.character(cut(x, 
                      cut_c2,
                      right=F,
                      labels=1:(length(cut_c2)-1))
  )
  
  x1
}



######################################## regroup_quartile ##################################
regroup_quartile=function(indata,var){
  value=unlist(indata[var])
  value_new=ifelse(value<quantile(value,na.rm = T)[2],1,
                   ifelse(value<quantile(value,na.rm = T)[3],2,
                          ifelse(value<quantile(value,na.rm = T)[4],3,4)))
  value_new
}

########################################  Test missing ########################################
output_nomissing_data=function(indata,var_all){
  indata1=indata[var_all]
  indata1$missing_index=apply(apply(indata1,2,function(x){ifelse(is.na(x),1,0)}),1,sum)
  outdata=indata[indata1$missing_index==0,]
  outdata
}

#x=complete.cases(indata[feature_in])
#indata=indata[x,]

######################################## cal dataset missing rate for dataset ##################################


cal_missing=function(indata,var_name_one){
  n_total=nrow(indata)
  value=unlist(indata[var_name_one])
  n_missing=n_total-table(is.na(value))[1]
  missing_rate=signif(round(n_missing/n_total,4),2)*100
  data.frame(
    var_name=var_name_one,
    missing_n_rate=paste0(n_missing, " (",missing_rate,"%)"),
    stringsAsFactors = F)
}

cal_missing_all=function(indata,var_name_all){
  bind_rows(lapply(var_name_all, function(var_name_one){
    cal_missing(indata,var_name_one)
  }
  ))
}
######################################## Exclude missing ##################################

cal_missing_ratio=function(x){
  round((length(x)-table(is.na(x))[1])/length(x),3)}


exclude_missing=function(indata,missing_cutoff){
  comp_need=names(indata[,-1])[apply(indata[,-1],2,cal_missing_ratio)<missing_cutoff]
  outdata=indata[c("studyid",comp_need)]
  outdata
}

######################################## Transformation ##################################
cal_log_natural=function(indata){
  log_natural_value=as.data.frame(apply(indata[,-1],2,log))
  log_natural_value1=cbind(studyid=indata[,1],log_natural_value)
  log_natural_value1
}


log_transform=function(indata,var_id){
  data_value=indata[-match(var_id,names(indata))]
  
  data_value1=as.data.frame(apply(data_value, 2, function(x){
    x1=x
    if(any(x<0)){x1=x-min(x,na.rm=T)}
    x1
  }))
  
  data_value_log=as.data.frame(apply(data_value1, 2, log))
  data_value_log1=cbind(indata[var_id],data_value_log)
  names(data_value_log1)[1]=var_id
  data_value_log1
}


######################################## Generate random value ##################################

generate_value_integer=function(n_need,max_value,min_value,type,sd_ratio=1){
  if(type=="right"){
    set.seed(1000)
    x0=round(rnorm(n_need*4, mean=0, sd=1*sd_ratio),3)
    x1=x0[x0>mean(x0)]
    x2=(x1-min(x1))/(max(x1)-min(x1))
    x3=x2*max_value
    x3=x3[x3>min_value]
    
    set.seed(0)
    if(length(x3)>=n_need){x4=sample(floor(x3),n_need)}
    set.seed(0)
    if(length(x3)<n_need){x4=sample(floor(x3),n_need,replace = T)}
    #table(x4)
    #hist(x4)
  }
  
  if(type=="left"){
    set.seed(1000)
    x0=round(rnorm(n_need*4, mean=0, sd=1*sd_ratio),3)
    x1=x0[x0<mean(x0)]
    x2=(x1-min(x1))/(max(x1)-min(x1))
    x3=x2*(max_value+1)
    x3=x3[x3>min_value]
    
    x3=ceiling(x3)
    x3=x3[!x3==(max_value+1)]
    
    set.seed(0)
    if(length(x3)>=n_need){x4=sample(x3,n_need)}
    set.seed(0)
    if(length(x3)<n_need){x4=sample(floor(x3),n_need,replace = T)}
    
    #hist(x4)
    #table(x4)
  }
  
  x4
}






######################################## Normalization ##################################
normalization_=function(indata,method_name){
  obj_meta=list(data=indata[,-1])
  normalized_=specmine::normalize(dataset = obj_meta,method =  method_name)
  normalized_value=normalized_$data
  outdata=as.data.frame(cbind(studyid=indata[,1],normalized_value))
  outdata
}

######################################## Scaling ##################################

scaling_=function(indata,method_name){
  scaled_=specmine::scaling_samples(datamat = indata[,-1],method =  method_name)
  outdata=as.data.frame(cbind(studyid=indata[,1],scaled_))
  outdata
}

# scaling --------------------------
# scale to [min_target, max_target]
scale_target=function(value,min_target,max_target,min_value,max_value){
  ifelse(value==-Inf | value== Inf,NA,
         ((value-min_value)/(max_value-min_value))*(max_target-min_target)+min_target
  )
}

