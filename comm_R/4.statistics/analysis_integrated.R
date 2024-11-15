########################################  Descriptive analysis for continuous variables ########################################


#cal_descriptive_con(indata,var_x_list,var_group,digit_n)

cal_descriptive_con_mean=function(indata,var_x_list,var_group,digit_n){
  bind_rows(lapply(var_x_list, function(var_x_one){
    bind_cols(
      cal_MeanSd_MultipleLevel_multipleVar(indata,var_x_one,var_group,digit_n),
      cal_annova_p(indata,var_group,var_x_one)[-1],
      Cal_trend_P_oneVar(indata,var_x_one,var_group)[-1]
    )
  }))
}


########################################  Descriptive analysis for categorical variables ########################################
#cal_descriptive_cat(indata,var_x_list,var_group,digit_n)

cal_descriptive_cat=function(indata,var_x_list,var_group,digit_n){
  
  x1=cal_freq_NoGroup(indata,var_x_list,digit_n)
  x2=cal_freq_set_cat(indata,var_group,var_x_list)[-(3:4)]
  x3=cal_chisq(var_x_list,var_group,indata)
  x4=cal_CochranArmitageTest_p(indata,var_group,var_x_list)
  xx=x1 %>% left_join(x2) %>% left_join(x3) %>% left_join(x4)
  xx
}

########################################  Descriptive analysis comprise of con and cat ########################################

cal_descriptive=function(indata,var_group,var_in_con,var_in_cat,data_label,digit_n){
  
  #Mean
  out_mean1=cal_MeanSd_MultipleLevel_multipleVar(indata,var_in_con,var_group,digit_n)
  out_mean2=bind_rows(lapply(var_in_con, function(var_in_con_one){
    cal_annova_p(indata,var_group,var_in_con_one)
  }))
  
  out_mean=out_mean1 %>% left_join(out_mean2) %>% mutate(variable_levels=NA,P_raw=NULL) %>% 
    left_join(Cal_trend_P_MultipleVar(indata,var_in_con,var_group))
  out=rbind(out_mean)
  
  var_in_cat_one=var_in_cat[1]
  
  #Freq
  if(!is.na(var_in_cat)[1]){
    out_freq1=bind_rows(lapply(var_in_cat, function(var_in_cat_one){
      freq1=cal_freq_one(indata,var_group,var_in_cat_one,"freq1")
      freq1$N_used=freq1$N_used[2]
      var_in_level_one=data_label$variable_levels[data_label$Variable==var_in_cat_one]
      freq1=freq1[freq1$variable_levels %in% as.character(var_in_level_one),] %>% 
        left_join(cal_chisq(var_in_cat_one,var_group,indata)) %>% 
        left_join(cal_CochranArmitageTest_p(indata,var_group,var_in_cat_one))
      freq1    
    }))
    
    out_freq2=bind_rows(lapply(var_in_cat, function(var_in_cat_one){
      freq2=cal_freq_NoGroup(indata,var_in_cat_one,digit_n)
      freq2$N_used=freq2$N_used[2]
      var_in_level_one=data_label$variable_levels[data_label$Variable==var_in_cat_one]
      freq2=freq2[freq2$variable_levels %in% var_in_level_one,]
      freq2
    }))
    
    out_freq=out_freq2 %>% left_join(out_freq1) %>% mutate(P_raw=NULL,index=NULL,test_name=NULL)
    out=rbind(out_freq)
  }
  
  #Set out
  if(!is.na(var_in_cat)[1]){out=rbind(out_mean,out_freq)}
  
  names(out)
  names(data_label)
  
  data_label$obs1=1:nrow(data_label)
  
  out1=out %>% left_join(data_label) %>% arrange(obs1)
  
  out1$Variable=out1$var_need_label
  
  out1$obs1=NULL
  out1$var_group=NULL
  out1$var_need_label=NULL
  out1$obs=NULL
  out1$variable_levels=NULL
  out1
}

########################################  Descriptive analysis comprise of con and cat with label dataset ########################################
#indata=data_follow_with_score1
#var_in_con=var_con
#var_in_cat=var_cat
#var_group="sex"
#data_label=data_label_basic

library(reshape2)

cal_descriptive_new=function(indata,var_in_con,var_in_cat,var_group,data_label){
  #Mean
  if(!is.na(var_in_con[1])){
    out_mean1=cal_MeanSd_MultipleLevel_multipleVar(indata,var_in_con,var_group,1)
    out_mean2=bind_rows(lapply(var_in_con, function(var_in_con_one){
      cal_annova_p(indata,var_group,var_in_con_one)
    }))
    out_mean=out_mean1 %>% left_join(out_mean2) %>% mutate(variable_levels=NA,P_raw=NULL) 
  }
  
  #Freq
  if(!is.na(var_in_cat[1])){
    out_freq1=bind_rows(lapply(var_in_cat, function(var_in_cat_one){
      cal_freq_NoGroup(indata,var_in_cat_one,1)
    }))
    out_freq2=bind_rows(lapply(var_in_cat, function(var_in_cat_one){
      cal_freq_chisq(indata,var_group,var_in_cat_one,"freq1")
    }))
    out_freq2$N_used=NULL
    out_freq=out_freq1 %>% left_join(out_freq2) %>% mutate(P_raw=NULL,index=NULL)
  }
  
  #Set out
  if(!is.na(var_in_con[1])){out=rbind(out_mean)}
  if(!is.na(var_in_cat[1])){out=rbind(out_freq)}
  if(!is.na(var_in_con[1]) & !is.na(var_in_cat[1])){out=rbind(out_freq,out_mean)}
  
  #out manipulation
  out$order1=1:nrow(out)
  out$order2=match(out$Variable,data_label$Variable)
  out=out %>% arrange(order2,order1) %>% mutate(order1=NULL, order2=NULL)
  
  
  out=out %>% left_join(data_label)
  out$variable_levels=ifelse(out$variable_levels==out$Variable,NA,out$variable_levels)
  out$Variable=data_label$var_need_label[match(out$Variable,data_label$Variable)]
  out$variable_levels=ifelse(is.na(out$variable_levels),out$Variable,out$variable_levels)
  out$variable_levels=ifelse(is.na(out$cat_level_label),out$variable_levels,out$cat_level_label)
  out[c("obs","var_need_label","var_group","cat_level_label")]=NULL
  out
}


















