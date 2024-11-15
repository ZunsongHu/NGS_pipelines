
######################################## logistic, one or multiple outcome, one or multiple x, with or without adjustments ##################################
#adj=""

cal_or=function(indata,var_outcome_list,var_x_list,adj,digit_n){
  style_number=paste0("%.",digit_n,"f")
  out_logit=
    bind_rows(lapply(var_outcome_list, function(var_outcome_one){
      
      bind_rows(lapply(var_x_list, function(var_x_one){
        indata$outcome=factor(unlist(indata[var_outcome_one]))
        
        if(is.na(adj)[1]|adj[1]==""){fit_formula=formula(paste0("outcome"," ~ ",paste0(c(var_x_one),collapse = "+")))}
        if(!(is.na(adj)[1]|adj[1]=="")){fit_formula=formula(paste0("outcome"," ~ ",paste0(c(var_x_one,adj),collapse = "+")))}
        
        logit=glm(fit_formula,family = "binomial",data=indata)
        
        out_one=data.frame(
          Outcome=var_outcome_one,
          N_used=logit$df.null+1,
          Cov=paste0(adj,collapse = "+"),
          Var_x=var_x_one,
          x_level=substr(names(coef(logit)),nchar(var_x_one)+1,nchar(names(coef(logit)))),
          Beta=sprintf(style_number,coef(logit)),
          OR=sprintf(style_number,exp(coef(logit))),
          CI=paste0(sprintf(style_number,exp(confint.default(logit)[,1])),"-",
                    sprintf(style_number,exp(confint.default(logit)[,2]))),
          P_raw=signif(summary(logit)$coefficients[,4],digits = 2),
          P=convert_p(summary(logit)$coefficients[,4]),
          stringsAsFactors = F
        )
        
        out_one=out_one[substr(names(coef(logit)),1,nchar(var_x_one))==var_x_one,]
        
        out_one  
      }))
    }))
  out_logit
}

######################################## logistic possion RR, one or multiple outcome, one or multiple x, with or without adjustments ##################################
