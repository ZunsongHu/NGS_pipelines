
########################################  2 group t test based on mean and sd ########################################
#value1="-0.11±0.233"
#value2="0.36±0.47"
#n1=60
#n2=60

#calculate_t_p_with_MeanSd(value1,value2,n1,n2)


calculate_t_p_using_MeanSd=function(value1,value2,n1,n2){
  mean1=as.numeric(word(value1,1,sep="±"))
  sd1=as.numeric(word(value1,2,sep="±"))
  
  
  mean2=as.numeric(word(value2,1,sep="±"))
  sd2=as.numeric(word(value2,2,sep="±"))
  
  mean_dif=mean1-mean2
  sp_squre=((n1-1)*(sd1^2)+(n2-1)*(sd2^2))/(n1+n2-2)
  se_dif=sqrt(sp_squre/n1+sp_squre/n2)
  
  
  t=mean_dif/se_dif
  P=2*pt(abs(t), df=(n1+n2-2),lower=F)
  
  return(list(t=t,P=P))
}




########################################  Description median ########################################
cal_Median_NogroupVariable=function(indata,var_list,digit_n){
  number_style=paste0("%.",digit_n,"f")
  bind_rows(lapply(var_list, function(var_one){
    data.frame(
      Variable=var_one,
      N_used=length(unlist(indata[var_one])[!is.na(unlist(indata[var_one]))]),
      All=
        paste0(sprintf(number_style,median(unlist(indata[var_one]),na.rm = T))," (",
               sprintf(number_style,quantile(unlist(indata[var_one]),0.25,na.rm = T))," - ",
               sprintf(number_style,quantile(unlist(indata[var_one]),0.75,na.rm = T)),")"
        ),
      stringsAsFactors = F
    )
  }))
}

cal_Median_WithGroupVariable=function(indata,var_list,var_group,digit_n){
  
  number_style=paste0("%.",digit_n,"f")
  
  value_group=unlist(indata[var_group])
  
  bind_rows(lapply(var_list, function(var_one){
    value_var=unlist(indata[var_one])
    
    out_one=data.frame(matrix(
      c(
        var_one,
        sum(!is.na(value_var) & !is.na(value_group)),
        paste0(
          sprintf(number_style,unlist(aggregate(value_var,list(value_group),function(x){median(x,na.rm = T)})[2]))," (",
          sprintf(number_style,unlist(aggregate(value_var,list(value_group),function(x){quantile(x,0.25,na.rm = T)})[2]))," - ",
          sprintf(number_style,unlist(aggregate(value_var,list(value_group),function(x){quantile(x,0.75,na.rm = T)})[2])),")"
        )
      ),
      nrow=1),stringsAsFactors = F)
    names(out_one)=c("Variable","N_used",paste0(var_group,"_",levels(factor(unlist(indata[,var_group])))))
    x1=kruskal.test(value_var~value_group)
    
    out_one$P_raw=x1$p.value
    out_one$P=convert_p(x1$p.value)
    out_one
  }))
}

cal_descriptive_con_median=function(indata,var_list,var_group,digit_n){
  x1=cal_Median_NogroupVariable(indata,var_list,digit_n)
  x2=cal_Median_WithGroupVariable(indata,var_list,var_group,digit_n)
  N_used=x2[c("Variable","N_used")]
  x1$N_used=NULL
  x2$N_used=NULL
  out_one=x1 %>% left_join(N_used) %>% left_join(x2)
  out_one
}





########################################  Calculate Mean and Sd of multiple var with No group ########################################
#digit_n=1

#cal_MeanSd_NogroupVariable(indata,outcome_con,1)

cal_MeanSd_NogroupVariable=function(indata,var_name,digit_n){
  mean_value=sprintf(paste0("%.",digit_n,"f"),apply(indata[var_name],2,function(x){mean(x,na.rm = T)}))
  sd_value=sprintf(paste0("%.",digit_n,"f"),apply(indata[var_name],2,function(x){sd(x,na.rm = T)}))
  
  median_value=sprintf(paste0("%.",digit_n,"f"),apply(indata[var_name],2,function(x){median(x,na.rm = T)}))
  q1_value=sprintf(paste0("%.",digit_n,"f"),apply(indata[var_name],2,function(x){quantile(x,0.25,na.rm = T)}))
  q3_value=sprintf(paste0("%.",digit_n,"f"),apply(indata[var_name],2,function(x){quantile(x,0.75,na.rm = T)}))
  
  min_value=sprintf(paste0("%.",digit_n,"f"),apply(indata[var_name],2,function(x){min(x,na.rm = T)}))
  max_value=sprintf(paste0("%.",digit_n,"f"),apply(indata[var_name],2,function(x){max(x,na.rm = T)}))
  
  N=apply(indata[var_name],2,function(x){sum(!is.na(x))})
  
  out=data.frame(
    Variable=var_name,
    N=N,
    Mean_sd=paste0(mean_value, " ± ",sd_value),
    
    Median_Q1_Q3=paste0(median_value, " (",q1_value," - ",q3_value,")"),
    
    Min_max=paste0(min_value,"-",max_value),
    
    stringsAsFactors = F)
  out
}

########################################  calculate Mean and Sd of one var with one group with multiple levels ########################################
cal_MeanSd_MultipleLevel_oneVar=function(indata,var_one,var_group,digit_n){
  out=data.frame(matrix(
    c(var_one,
      sum(aggregate(indata[var_one],by=list(unlist(indata[var_group])),function(x){sum(!is.na(x))})[2]),
      paste0(
        sprintf(paste0("%.",digit_n,"f"),mean(unlist(indata[!is.na(indata[var_group]),var_one]),na.rm=T)),
        " ± ",
        sprintf(paste0("%.",digit_n,"f"),sd(unlist(indata[!is.na(indata[var_group]),var_one]),na.rm=T))
      ),
      paste0(
        t(aggregate(indata[var_one],by=list(unlist(indata[var_group])),function(x){sprintf(paste0("%.",digit_n,"f"),mean(x,na.rm = T))})[,-1]),
        " ± ",
        t(aggregate(indata[var_one],by=list(unlist(indata[var_group])),function(x){sprintf(paste0("%.",digit_n,"f"),sd(x,na.rm = T))})[,-1])
      )
    ),
    nrow=1
  ),
  stringsAsFactors = F)
  
  names(out)=c("Variable","N_used","All",paste0(var_group,"_",levels(factor(unlist(indata[,var_group])))))
  out$N_used=as.numeric(out$N)
  out
}

########################################  calculate Mean and Sd of multiple var with one group with multiple levels ########################################

cal_MeanSd_MultipleLevel_multipleVar=function(indata,var_multiple,var_group,digit_n){
  out1=bind_rows(
    lapply(var_multiple, function(var_one){
      cal_MeanSd_MultipleLevel_oneVar(indata,var_one,var_group,digit_n)
    }
    )
  )
  out1  
}



########################################  Calculate Mean and Sd of multiple var with one group variable and one strata variable ########################################

cal_MeanSd_oneGroup_oneStrata=function(indata,var_x,var_strata,var_group,digit_n){
  number_style=paste0("%.",digit_n,"f")
  value_strata_level=levels(factor(unlist(indata[var_strata])))
  
  cal_meansd_one_temp=function(i){
    indata1=indata[unlist(indata[var_strata])==value_strata_level[i],]
    N_used=length(unlist(indata1[var_x])[!is.na(unlist(indata1[var_x]))])
    
    mean_all=mean(unlist(indata1[var_x]),na.rm=T)
    sd_all=sd(unlist(indata1[var_x]),na.rm=T)
    All=paste0(sprintf(number_style,mean_all)," ± ",sprintf(number_style,sd_all))
    
    mean_g=aggregate(indata1[var_x],by=list(unlist(indata1[var_group])),function(x){mean(x,na.rm=T)})
    sd_g=aggregate(indata1[var_x],by=list(unlist(indata1[var_group])),function(x){sd(x,na.rm=T)})
    mean_level=paste0(sprintf(number_style,mean_g[,2])," ± ",sprintf(number_style,sd_g[,2]))
    
    out_one=data.frame(c(
      var_strata,
      value_strata_level[i],
      N_used,All,mean_level,
      cal_annova_p_one(indata1,var_group,var_x)[c(2,3)],
      Cal_trend_P_oneVar(indata1,var_x,var_group)[c(2,3)]),
      stringsAsFactors = F
    )
    
    names(out_one)=c("Strat",'Strat_levels','N_used','All',
                     mean_g[,1],
                     'P_raw','P','P_trend_raw','P_trend'
    )
    out_one
  }
  
  label=data.frame(Strat=var_strata,stringsAsFactors = F)
  out_mean=bind_rows(label,bind_rows(lapply(1:length(value_strata_level), cal_meansd_one_temp)))
  out_mean
}

























