


#' cal_freqPer_2groups
#'
#' @param df_in 
#' @param level_high 
#' @param level_low 
#'
#' @return
#' @export
#'
#' @examples

#Sum of percentage in each level of big group = 1
#var_passing must be sub level(s) of level_high
cal_freqPer_2groups=function(df_in,level_low,level_high,var_passing=NULL){
  df_in1=df_in[c(level_high,level_low)]
  names(df_in1)=c('level_high','level_low')
  
  # df_in1$level_high=as.character(df_in1$level_high)
  # df_in1$level_low=as.character(df_in1$level_low)
  
  df_in2=table(df_in1$level_low,df_in1$level_high) %>% as.data.frame() %>% 
    dplyr::rename(level_low = Var1, level_high=Var2,n_g_level = Freq) %>% 
    left_join(table(df_in1$level_high) %>% as.data.frame() %>% dplyr::rename(level_high = Var1, n_g = Freq)) %>% 
    mutate(per=n_g_level/n_g,per1=round(per*100,1),per2=paste0(per1,"%")) %>% 
    distinct()  
  
  names(df_in2)[1:2]=c(c(level_low,level_high))
  
  if(!is.null(var_passing)){
    df_passing=df_in[c(level_high,var_passing)] %>% distinct()
    df_in2=df_in2 %>% left_join(df_passing)
  }
  
  df_in2
}

#anchor=mid: Sum of percentage in each level of big group and mid group = 1
#anchor=high: Sum of percentage in each level of big group = 1
#anchor=mid: Sum of percentage in each level of mid group = 1

cal_freqPer_3groups=function(df_in,level_high,level_mid,level_low,anchor="mid"){
  if(!all(c(level_low,level_mid,level_high) %in% names(df_in))){stop("Var(s) not in dataset")}
  
  df_in1=df_in[c(level_low,level_mid,level_high)]
  names(df_in1)=c('level_low',"level_mid",'level_high')
  
  if(anchor=="mid"){
    df_in2=df_in1 %>% 
      group_by(level_high,level_mid) %>% mutate(n_g=n()) %>% 
      group_by(level_high,level_mid,level_low) %>% 
      mutate(n_g_level=n(),
             per=round(n_g_level/n_g,6)) %>% 
      select(level_low,level_mid,level_high,n_g,n_g_level,per) %>% distinct()
  }
  
  if(anchor=="high"){
    df_in2=df_in1 %>% 
      group_by(level_high) %>% mutate(n_g=n()) %>% 
      group_by(level_high,level_mid,level_low) %>% 
      mutate(n_g_level=n(),
             per=round(n_g_level/n_g,6)) %>% 
      select(level_low,level_mid,level_high,n_g,n_g_level,per) %>% distinct()
  }
  
  names(df_in2)[1:3]=c(level_low,level_mid,level_high)
  df_in2
  
}

#calculate freqs --------------
cal_freq_NoGroup=function(indata,var_list,digit_n){
  bind_rows(lapply(var_list, function(var_x){
    value_x=unlist(indata[var_x])
    out1=data.frame(table(value_x))
    out2=data.frame(
      Variable=var_x,
      variable_levels=paste0(as.character(unlist(out1[,1]))),
      N_used=as.character(sum(!is.na(value_x))),
      All=
        paste0(unlist(out1[,2]),
               " (",
               sprintf("%.1f",unlist(out1[,2])/sum(unlist(out1[,2]))*100),
               "%)"
        ),
      freq=unlist(out1[,2]),
      percentage=paste0(sprintf("%.1f",unlist(out1[,2])/sum(unlist(out1[,2]))*100),"%"),
      stringsAsFactors = F
    )
    out2$N_used[2:nrow(out2)]=NA
    bind_rows(data.frame(Variable=var_x,variable_levels=var_x,stringsAsFactors = F),out2) 
  }))
  
}




########################################  Calculate freq, with group variable, one variable ########################################

#cal_freq_one(indata,group,var_cat_one,"freq1")

#group="all"

# calculation for one group variable
cal_freq_one=function(indata,group,var_cat_one,freq_need="freq1"){
  value_group=unlist(indata[group])
  value_cat_one=unlist(indata[var_cat_one])
  
  factor_group=factor(value_group)
  factor_cat=factor(value_cat_one)
  
  levels(factor_group)[length(levels(factor_group))]
  levels(factor_cat)[length(levels(factor_cat))]
  
  freq_raw=table(value_cat_one,value_group)
  freq_sum_col=apply(freq_raw,2,sum)
  
  #0
  freq1=data.frame(Variable=rep(var_cat_one,length(levels(factor_cat))+1),
                   variable_levels=c(var_cat_one,row.names(freq_raw)),
                   index=c(NA,paste0(var_cat_one,"_",row.names(freq_raw))),
                   stringsAsFactors = F)
  freq1$N_used=NA;freq1$N_used[2]=sum(freq_raw)
  
  #1
  if (freq_need=="freq1"){
    percentage_raw1=freq_raw
    for (i in 1:length(levels(factor_cat))){
      percentage_raw1[i,]=percentage_raw1[i,]/freq_sum_col
    }
    cal=function(x){paste0(sprintf("%3.3f",round(x,10)*100),"%")}
    percentage_raw1[]=cal(percentage_raw1)
    percentage_raw1[]=paste0(freq_raw," (",percentage_raw1,")")
    percentage_raw1=as.data.frame(percentage_raw1,stringsAsFactors = F)
    percentage_raw1=dcast(percentage_raw1,value_cat_one~value_group,value.var = 'Freq')
    
    names(percentage_raw1)=c("variable_levels",paste0(group,"_",levels(factor_group)))
    freq_out=freq1 %>% left_join(percentage_raw1)
    
  }
  #2
  if (freq_need=="freq2"){
    cal=function(x){paste0(sprintf("%3.3f",round(x,10)*100),"%")}
    percentage_raw2=freq_raw/sum(freq_raw)
    percentage_raw2[]=cal(percentage_raw2)
    percentage_raw2[]=paste0(freq_raw," (",percentage_raw2,")")
    percentage_raw2=as.data.frame(percentage_raw2,stringsAsFactors = F)
    percentage_raw2=dcast(percentage_raw2,value_cat_one~value_group,value.var = 'Freq')
    names(percentage_raw2)=c("variable_levels",paste0(group,"_",levels(factor_group)))
    freq_out=freq1 %>% left_join(percentage_raw2)
  }
  freq_out[]=lapply(freq_out, as.character)
  freq_out
}


#freq_base1=cal_freq_set_cat(data_merge,'all',var_cat_base)
cal_freq_set_cat=function(indata,group,var_cat_all){
  freq_out_set=
    bind_rows(lapply(var_cat_all,function(var_cat_one){
      cal_freq_one(indata,group,var_cat_one,"freq1")
    }))
  freq_out_set
}

########################################  Calculate chisq ########################################
#indata=data_1
#var_group="V3"
#var_freq="V1"
cal_chisq=function(var_freq_list,var_group,indata){
  freq_group  = factor(unlist(indata[var_group]))
  
  bind_rows(lapply(var_freq_list, function(var_freq){
    print(var_freq)
    freq_var    = factor(unlist(indata[var_freq]))
    
    n_min=min(table(freq_group,freq_var))
    
    x1=data.frame(table(freq_group,freq_var))
    x1$value=paste0("value_",x1$freq_var)
    x2=tidyr::spread(x1[-2], value, Freq)
    
    if(any(apply(x2[2:ncol(x2)],1,sum)==0) |  any(apply(x2[2:ncol(x2)],2,sum)==0)){
      chisq_p=NA
    } else{
      
      if(n_min>4){
        chisq_p=chisq.test(freq_var,freq_group)$p.value}
      
      if(n_min<=4){
        chisq_p=fisher.test(table(freq_group,freq_var),simulate.p.value = T,alternative = "two.sided")$p.value
      }
    }
    data.frame(
      Variable=var_freq,
      P_raw=signif(chisq_p,digits = 6),
      P=convert_p(chisq_p),
      stringsAsFactors = F)
  }))
  
}





#calculate freqs--------------
cal_freq_NoGroup=function(indata,var_list,digit_n){
  bind_rows(lapply(var_list, function(var_x){
    value_x=unlist(indata[var_x])
    out1=data.frame(table(value_x))
    out2=data.frame(
      Variable=var_x,
      variable_levels=paste0(as.character(unlist(out1[,1]))),
      N_used=as.character(sum(!is.na(value_x))),
      All=
        paste0(unlist(out1[,2]),
               " (",
               sprintf("%.1f",unlist(out1[,2])/sum(unlist(out1[,2]))*100),
               "%)"
        ),
      freq=unlist(out1[,2]),
      percentage=paste0(sprintf("%.1f",unlist(out1[,2])/sum(unlist(out1[,2]))*100),"%"),
      stringsAsFactors = F
    )
    out2$N_used[2:nrow(out2)]=NA
    bind_rows(data.frame(Variable=var_x,variable_levels=var_x,stringsAsFactors = F),out2) 
  }))
  
}

########################################  Calculate freq, no group variable, one or multiple variables ########################################

#value_x=unlist(indata["hp"])
#var_x="hp"

cal_freq_NoGroup=function(indata,var_list,digit_n){
  bind_rows(lapply(var_list, function(var_x){
    value_x=unlist(indata[var_x])
    out1=data.frame(table(value_x))
    out2=data.frame(
      Variable=var_x,
      variable_levels=paste0(as.character(unlist(out1[,1]))),
      N_used=as.character(sum(!is.na(value_x))),
      All=
        paste0(unlist(out1[,2]),
               " (",
               sprintf("%.1f",unlist(out1[,2])/sum(unlist(out1[,2]))*100),
               "%)"
        ),
      stringsAsFactors = F
    )
    out2$N_used[2:nrow(out2)]=NA
    bind_rows(data.frame(Variable=var_x,variable_levels=var_x,stringsAsFactors = F),out2) 
  }))
  
}


########################################  Calculate freq, with group variable, one variable ########################################

#############calculation for one group variable
cal_freq_one=function(indata,group,var_cat_one,freq_need="freq1"){
  value_group=unlist(indata[group])
  value_cat_one=unlist(indata[var_cat_one])
  
  factor_group=factor(value_group)
  factor_cat=factor(value_cat_one)
  
  levels(factor_group)[length(levels(factor_group))]
  levels(factor_cat)[length(levels(factor_cat))]
  
  freq_raw=table(value_cat_one,value_group)
  freq_sum_col=apply(freq_raw,2,sum)
  
  #0
  freq1=data.frame(Variable=rep(var_cat_one,length(levels(factor_cat))+1),
                   variable_levels=c(var_cat_one,row.names(freq_raw)),
                   index=c(NA,paste0(var_cat_one,"_",row.names(freq_raw))),
                   stringsAsFactors = F)
  freq1$N_used=NA;freq1$N_used[2]=sum(freq_raw)
  
  #1
  if (freq_need=="freq1"){
    percentage_raw1=freq_raw
    for (i in 1:length(levels(factor_cat))){
      percentage_raw1[i,]=percentage_raw1[i,]/freq_sum_col
    }
    cal=function(x){paste0(sprintf("%3.3f",round(x,10)*100),"%")}
    percentage_raw1[]=cal(percentage_raw1)
    percentage_raw1[]=paste0(freq_raw," (",percentage_raw1,")")
    percentage_raw1=as.data.frame(percentage_raw1,stringsAsFactors = F)
    percentage_raw1=dcast(percentage_raw1,value_cat_one~value_group,value.var = 'Freq')
    
    names(percentage_raw1)=c("variable_levels",paste0(group,"_",levels(factor_group)))
    freq_out=freq1 %>% left_join(percentage_raw1)
    
  }
  #2
  if (freq_need=="freq2"){
    cal=function(x){paste0(sprintf("%3.3f",round(x,10)*100),"%")}
    percentage_raw2=freq_raw/sum(freq_raw)
    percentage_raw2[]=cal(percentage_raw2)
    percentage_raw2[]=paste0(freq_raw," (",percentage_raw2,")")
    percentage_raw2=as.data.frame(percentage_raw2,stringsAsFactors = F)
    percentage_raw2=dcast(percentage_raw2,value_cat_one~value_group,value.var = 'Freq')
    names(percentage_raw2)=c("variable_levels",paste0(group,"_",levels(factor_group)))
    freq_out=freq1 %>% left_join(percentage_raw2)
  }
  freq_out[]=lapply(freq_out, as.character)
  freq_out
}


#freq_base1=cal_freq_set_cat(data_merge,'all',var_cat_base)
cal_freq_set_cat=function(indata,group,var_cat_all){
  freq_out_set=
    bind_rows(lapply(var_cat_all,function(var_cat_one){
      cal_freq_one(indata,group,var_cat_one,"freq1")
    }))
  freq_out_set
}


########################################  Calculate chisq ########################################
#indata=data_1
#var_group="V3"
#var_freq="V1"
cal_chisq=function(var_freq_list,var_group,indata){
  freq_group  = factor(unlist(indata[var_group]))
  
  bind_rows(lapply(var_freq_list, function(var_freq){
    print(var_freq)
    freq_var    = factor(unlist(indata[var_freq]))
    
    n_min=min(table(freq_group,freq_var))
    
    x1=data.frame(table(freq_group,freq_var))
    x1$value=paste0("value_",x1$freq_var)
    x2=tidyr::spread(x1[-2], value, Freq)
    
    if(any(apply(x2[2:ncol(x2)],1,sum)==0) |  any(apply(x2[2:ncol(x2)],2,sum)==0)){
      chisq_p=NA
    } else{
      
      if(n_min>4){
        chisq_p=chisq.test(freq_var,freq_group)$p.value}
      
      if(n_min<=4){
        chisq_p=fisher.test(table(freq_group,freq_var),simulate.p.value = T,alternative = "two.sided")$p.value
      }
    }
    data.frame(
      Variable=var_freq,
      P_raw=signif(chisq_p,digits = 6),
      P=convert_p(chisq_p),
      stringsAsFactors = F)
  }))
  
}

######################################## Merge freq and chisq, one or multiple x var ########################################

#indata=indata
#group="source1"
#var_cat_one="points_sum"

#cal_freq_chisq_one(indata,"source1","points_sum","freq1")


cal_freq_chisq=function(indata,group,var_cat_list,freq_type){
  bind_rows(lapply(var_cat_list, function(var_cat_one){
    out_freq=cal_freq_one(indata,group,var_cat_one,freq_type)
    out_chisq=cal_chisq(var_cat_one,group,indata)
    out_freq=out_freq %>% left_join(out_chisq)
    n=1:nrow(out_freq)
    out_freq$P=ifelse(n==2,out_freq$P,NA)
    out_freq
  }))
}


########################################  Extract Percentage ########################################

extract_percentage=function(x){
  perc=grepl("%",x)
  x[perc]=word(word(x[perc],2,sep="[(]"),1,sep="[)]")
  x
}

########################################  Extract freq ########################################
extract_freq=function(x){
  as.numeric(word(x,1,sep=" "))
}




















