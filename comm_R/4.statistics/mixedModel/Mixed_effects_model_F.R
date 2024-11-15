######################################## linear mixed effects model ##################################
#indata=data_avidity

#names(indata)

#outcome="value"
#co_var=c("Group","time")
#random_var=c("ID")

cal_mixed_model_one=function(indata,outcome,co_var,random_var){
  
  fit_formula=formula(
    paste0(outcome," ~  ",paste0(co_var,collapse = "+"),"+",paste0(paste0("(1|",random_var,")"),collapse = "+"))
  )
  
  fit_mixed <- lmerTest::lmer(fit_formula,data =  indata,REML=F)
  
  coef_fixed=as.data.frame(summary(fit_mixed)$coefficients)
  
  cbind(
    data.frame( Variables=row.names(coef_fixed),stringsAsFactors = F),
    coef_fixed
  )
  
}

#Mixed_out_avidity=cal_mixed_model_one(data_avidity,"value",c("Group","time"),c("ID"))

######################################## cox mixed effects model ##################################


#install.packages("coxme")



cal_cox_mixed=function(indata,outcome,outcome_t,var_id,var_x,var_adj,digit_n){
  library(survival)
  style_number=paste0("%.",digit_n,"f")
  
  if(!is.na(var_adj)[1]){var_x_in=c(var_adj,var_x)}
  if(is.na(var_adj)[1]){var_x_in=c(var_x)}
  random_items=paste0("(1|",var_id,")")
  
  formula_mixed=formula(paste0(paste0("Surv(",outcome_t,",",outcome,")"),"~",paste0(c(var_x_in,random_items),collapse = "+")))
  
  cox_mixed=coxme::coxme(formula_mixed, data = indata)
  #cox_mixed
  sink("mixed_temp_out.txt")
  print(cox_mixed)
  sink()
  
  cox_mixed1=readLines("mixed_temp_out.txt")
  cox_mixed2=cox_mixed1[((1:length(cox_mixed1))[cox_mixed1=="Fixed coefficients"]+2):
                          ((1:length(cox_mixed1))[cox_mixed1=="Random effects"]-2) ]
  file.remove("mixed_temp_out.txt")
  
  out_one=data.frame(
    #variable=var_x_in,
    variable_levels=word(cox_mixed2,1,sep="[ ]+"),
    beta=as.numeric(cox_mixed$coefficients),
    se=as.numeric(word(cox_mixed2,4,sep="[ ]+")),
    p_raw=as.numeric(word(cox_mixed2,6,sep="[ ]+")),
    stringsAsFactors = F
  )
  
  out_one1=out_one %>% mutate(RR=sprintf(style_number,exp(beta)),
                              rr_ll=sprintf(style_number,exp(beta-1.96*se)),
                              rr_ul=sprintf(style_number,exp(beta+1.96*se))
  ) %>% transmute(
    outcome=outcome,
    #variable=variable,
    variable_levels=variable_levels,
    #levels=substr(variable_levels,(nchar(variable)+1),nchar(variable_levels)),
    N_used=cox_mixed$n[2],
    beta=sprintf(style_number,beta),
    se=sprintf(style_number,se),
    RR=paste0(RR," (",rr_ll," - ",rr_ul,")"),
    P_raw=p_raw
  )
  return(out_one1)
}

















