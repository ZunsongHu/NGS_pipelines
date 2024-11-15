########################################  calculate_linear_fitStatistic ########################################

#indata=data_nutrition
#outcome="grni"
#var_in=unlist(data_list_one[1,])

calculate_linear_fitStatistic=function(indata,outcome,var_in){
  indata1=na.omit(indata[c(outcome,var_in)])
  value_outcome=unlist(indata1[outcome])
  fit_formula=formula(paste0(outcome,"~",paste0(c(var_in),collapse = "+")))
  
  fit1=summary(glm(fit_formula,data=indata1))
  fit2=summary(lm(fit_formula,data=indata1))
  value_predict=predict(glm(fit_formula,data=indata1))
  
  data_one=data.frame(
    aic=fit1$aic,
    r_square=fit2$r.squared,
    r_square_adj=fit2$adj.r.squared,
    P_Ftest=pf(fit2$fstatistic[1] , (fit2$df[1]-1), fit2$df[2], lower.tail = FALSE),
    mae=Metrics::mae(value_outcome, value_predict),
    rmse=Metrics::rmse(value_outcome, value_predict),
    corr_original_predicted=cor(value_outcome,value_predict),
    stringsAsFactors = F
  )
  data_one
}




########################################  Convert p values ########################################
convert_p=function(P_raw){
  P=ifelse(P_raw<0.001,"<0.001",
           ifelse(P_raw<0.01,sprintf("%.3f",P_raw),
                  ifelse(P_raw<0.045 | P_raw==0.045,sprintf("%.2f",P_raw),
                         ifelse(P_raw<0.055 | P_raw==0.055,sprintf("%.3f",P_raw),
                                ifelse(P_raw>0.055,sprintf("%.2f",P_raw),NA
                                )))))
  P
}


######################################## Paired compare for each group against the lowest group  ##################################
cal_pairedCompare_P_oneVar=function(indata,var_one,var_group){
  value_levels=levels(factor(unlist(indata[var_group])))
  value_levels_base=value_levels[1]
  value_levels_comp=value_levels[-1]
  value_levels_comp_one=value_levels_comp[1]
  
  out1=bind_cols(
    lapply(value_levels_comp, function(value_levels_comp_one){
      indata1=indata[unlist(indata[var_group]) %in% c(value_levels_comp_one,value_levels_base),]
      
      fit_one=summary(lm(unlist(indata1[var_one])~unlist(indata1[var_group])))
      P_raw=fit_one$coefficients[2,4]
      
      out1=data.frame(
        
        P=ifelse(P_raw<0.001,"<0.001",
                 ifelse(P_raw<0.01,sprintf("%.3f",P_raw),
                        ifelse(P_raw<0.045 | P_raw==0.045,sprintf("%.2f",P_raw),
                               ifelse(P_raw<0.055 | P_raw==0.055,sprintf("%.3f",P_raw),
                                      ifelse(P_raw>0.055,sprintf("%.2f",P_raw),NA
                                      ))))),
        stringsAsFactors = F
      )
      names(out1)=paste0("P ",paste0(c(value_levels_comp_one,value_levels_base),collapse = " VS "))
      out1
      
    }
    )
  )
  out2=cbind(Variable=var_one,out1)
  out2$Variable=as.character(out2$Variable)
  out2
}


cal_pairedCompare_P_MultipleVar=function(indata,var_multiple,var_group){
  out1=bind_rows(
    lapply(var_multiple, function(var_one){
      cal_pairedCompare_P_oneVar(indata,var_one,var_group)
    })
  )
  out1
}

######################################## Trend test based on median of each group  ##################################
#cal_MeanSd_MultipleLevel_multipleVar(data_folate,outcome_con,'M1PFolate_result_4g',number_digital)

Cal_trend_P_oneVar=function(indata,var_one,var_group){
  #print(var_one)
  value_levels=levels(factor(unlist(indata[var_group])))
  value_median=aggregate(indata[var_one],by = list(unlist(indata[var_group])),function(x){median(x,na.rm=T)})[,2]
  value_for_fit=value_median[match(unlist(indata[var_group]),value_levels)]
  order_for_fit=c(1:length(value_levels))[match(unlist(indata[var_group]),value_levels)]
  fit_trend=summary(lm(unlist(indata[var_one])~value_for_fit))
  
  P=NA
  
  if(nrow(fit_trend$coefficients)>=2){
    P_raw=fit_trend$coefficients[2,4]
    
    P=ifelse(P_raw<0.001,"<0.001",
             ifelse(P_raw<0.01,sprintf("%.3f",P_raw),
                    ifelse(P_raw<0.045 | P_raw==0.045,sprintf("%.2f",P_raw),
                           ifelse(P_raw<0.055 | P_raw==0.055,sprintf("%.3f",P_raw),
                                  ifelse(P_raw>0.055,sprintf("%.2f",P_raw),NA
                                  )))))
  }
  
  if(nrow(fit_trend$coefficients)==1){
    fit_trend=summary(lm(unlist(indata[var_one])~order_for_fit))
    P_raw=fit_trend$coefficients[2,4]
    P=ifelse(P_raw<0.001,"<0.001",
             ifelse(P_raw<0.01,sprintf("%.3f",P_raw),
                    ifelse(P_raw<0.045 | P_raw==0.045,sprintf("%.2f",P_raw),
                           ifelse(P_raw<0.055 | P_raw==0.055,sprintf("%.3f",P_raw),
                                  ifelse(P_raw>0.055,sprintf("%.2f",P_raw),NA
                                  )))))
  }
  
  
  data.frame(
    Variable=var_one,
    P_trend_raw=signif(P_raw,2),
    P_trend=P,
    stringsAsFactors = F
  )
}


Cal_trend_P_MultipleVar=function(indata,var_multiple,var_group){
  out1=bind_rows(
    lapply(var_multiple, function(var_one){
      Cal_trend_P_oneVar(indata,var_one,var_group)
    }
    )
  )
  out1
}

######################################## Partial Correlation analysis  ##################################
#install.packages("ppcor")
#library(ppcor)

#outcome_one=var_for_analysis[1]
#adj=var_adj
#var_one=var_diversity[1]
#indata=Merge_diversity
#library(ppcor)

#adj=""


#install.packages("ppcor")

cal_pcor_one=function(indata,outcome_one,var_one,adj){
  if(!adj[1]==""){var_all=c(outcome_one,adj,var_one)}
  if(adj[1]==""){var_all=c(outcome_one,var_one)}
  indata1=indata[var_all]
  indata2=indata1[!apply(indata1,1,function(x){any(is.na(x))}),]
  
  if(!adj[1]==""){
    pcor_fit=ppcor::pcor.test(x = indata2[,var_one], y = indata2[,outcome_one], z = indata2[,adj],method = "spearman")
    out=data.frame(
      Outcome=outcome_one,
      Variable=var_one,
      N=pcor_fit$n,
      Estimate=round(pcor_fit$estimate,5),
      P=signif(pcor_fit$p.value,2),
      Method=as.character(pcor_fit$Method),
      Adjustment=paste0(adj,collapse = "+"),
      stringsAsFactors = F)
  }
  if(adj[1]==""){
    cor_fit=cor.test(unlist(indata[outcome_one]),unlist(indata[var_one]),method = "spearman",exact = F)
    out=data.frame(
        Outcome=outcome_one,
        Variable=var_one,
        N=NA,
        Estimate=round(cor_fit$estimate,5),
        P=signif(cor_fit$p.value,2),
        Method=as.character(cor_fit$method),
        Adjustment="",
        stringsAsFactors = F)
    
  }
  out
}

cal_pcor_set=function(indata,outcome_one,var_list,adj){
  out_pcor=bind_rows(lapply(var_list, function(var_one){
    cal_pcor_one(indata,outcome_one,var_one,adj)
  }))
  out_pcor
}



######################################## Regular linear regression, one outcome or outcome list, one x or x list, with or without adjustments ##################################
linear_regression=function(indata,outcome_list,var_x_list,adj,digit_n){
  
  out_one=
    bind_rows(lapply(outcome_list, function(outcome_one){
      
      bind_rows(lapply(var_x_list, function(var_one){
        print(paste0(outcome_one,"~",var_one))
        
        fit_formula=formula(paste0(outcome_one,"~",paste0(c(var_one,adj),collapse = "+")))
        
        if(is.na(adj)[1]|adj[1]==""){fit_formula=formula(paste0(outcome_one," ~ ",paste0(c(var_one),collapse = "+")))}
        if(!(is.na(adj)[1]|adj[1]=="")){fit_formula=formula(paste0(outcome_one," ~ ",paste0(c(var_one,adj),collapse = "+")))}
        
        fit_lm=summary(lm(fit_formula,data=indata))
        
        value_coefficients=as.data.frame(fit_lm$coefficients)
        value_coefficients=value_coefficients[substr(row.names(value_coefficients),1,nchar(var_one))==var_one,]
        
        out_one=data.frame(
          outcome=outcome_one,
          N_used=sum(fit_lm$df[1:2]),
          cov=paste0(adj,collapse = "+"),
          
          var_x=var_one,
          x_level=substr(row.names(value_coefficients),nchar(var_one)+1,nchar(row.names(value_coefficients))),
          
          Beta_se=paste0(
            sprintf(paste0("%.",digit_n,"f"),round(value_coefficients[,1],4))," (",
            sprintf(paste0("%.",digit_n,"f"),round(value_coefficients[,2],4)),")"),
          standard_beta=value_coefficients[,1]/value_coefficients[,2],
          P_raw=value_coefficients[,4],
          P=convert_p(value_coefficients[,4]),
          stringsAsFactors = F
        )
        
        out_one
      }))
    }))
  out_one
}

######################################## calculate multiple, out all beta ##################################
calculate_mulpitle_betaSe=function(indata,outcome_one,var_in,digit_n){
  indata1=na.omit(indata[c(outcome_one,var_in)])
  
  digit_style=paste0("%.",digit_n,"f")
  
  fit_formula=formula(paste0(outcome_one," ~ ",paste0(var_in,collapse = "+")))
  
  fit_lm=glm(fit_formula,data=indata1)
  fit_summary=summary(fit_lm)
  
  coef=fit_summary$coefficients
  out_one=
    data.frame(
      outcome=outcome_one,
      var_x=row.names(coef),
      Beta=sprintf(digit_style,coef[,1]),
      se=sprintf(digit_style,coef[,2]),
      standard_beta=coef[,1]/coef[,2],
      p_raw=coef[,4],
      P=convert_p(coef[,4]),
      stringsAsFactors = F
    )
  out_one[!out_one$var_x=="(Intercept)",]
}


######################################## calculate predicted ##################################



calculate_predict=function(indata,outcome_one,var_in,var_id){
  indata1=na.omit(indata[c(var_id,outcome_one,var_in)])
  
  fit_formula=formula(paste0(outcome_one," ~ ",paste0(var_in,collapse = "+")))
  
  fit_lm=glm(fit_formula,data=indata1)
  fit_summary=summary(fit_lm)
  out_one=data.frame(var_id=unlist(indata1[var_id]),
                     predicted=predict(fit_lm),
                     stringsAsFactors = F
  )
  names(out_one)=c(var_id,paste0(outcome_one,"_predicted"))
  out_one
  
}



######################################## Regular correlation one ##################################
#indata=data_pheno_need_MI450 %>% left_join(out_knn_median_auto)
#outcome=feature_in_all[1]
#var_x_one="Energy_intake_adj"

cal_cor_one=function(indata,outcome,var_x_one){
  out_one=data.frame(
    outcome=outcome,
    var_x=var_x_one,
    cor_pearson=cor(unlist(indata[outcome]),unlist(indata[var_x_one]),method = "pearson",use = "complete.obs"),
    cor_spearman=cor(unlist(indata[outcome]),unlist(indata[var_x_one]),method = "spearman",use = "complete.obs"),
    P_spearman=signif(cor.test(unlist(indata[outcome]),unlist(indata[var_x_one]),method = "spearman",exact = F)$p.value,2),
    stringsAsFactors = F
  )
  out_one
}

#set when outcome is more than x numbers
#outcome_all=feature_in_all
#var_x_all=c('mother_age',"Pre_BMI",'FastFodd_pattern','Energy_intake_adj','Protein_intake','Fat_intake','Carbon_intake')

cal_cor_set=function(indata,outcome_all,var_x_all){
  cor_out=bind_cols(
    lapply(var_x_all, function(var_x_one){
      print(var_x_one)
      cor_set=bind_rows(
        lapply(outcome_all, function(outcome){
          cal_cor_one(indata,outcome,var_x_one)
        }
        ))
      names(cor_set)=paste0(names(cor_set),"_",var_x_one)
      cor_set[c(3,4)]
    }
    ))
  out=cbind(outcome_all,cor_out)
  out$outcome_all=as.character(out$outcome_all)
  out
}

#cor_out=cal_cor_set(indata,outcome_all,var_x_all)






