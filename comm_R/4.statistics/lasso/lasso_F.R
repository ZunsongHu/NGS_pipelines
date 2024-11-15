#data_beta should have tow columns, the name of variables(Variables) and the beta values(Beta).

cal_score_basedOnBeta=function(indata,data_beta){
  score=0
  for (j in 1:length(data_beta$Variables)){
    score=score+unlist(indata[data_beta$Variables[j]])*(data_beta$Beta[j])}
  score
}

########################################### lasso cox, calculate lambda, beta and predicted value ##################################################


lasso_cox_one=function(indata,outcome,outcome_t,var_x,seed){
  fit_formula=formula(paste0(" ~ ",paste(var_x,collapse = "+")))
  x=model.matrix(fit_formula, indata)
  y=Surv(unlist(indata[outcome_t]), unlist(indata[outcome]))
  
  #fit <- glmnet(x, y, family="cox", alpha=1)
  set.seed(seed)
  cvfit = glmnet::cv.glmnet(x, y,family="cox", alpha=1,nfolds = 10)
  #plot(cvfit)
  b_fit_cv=coef(cvfit,s="lambda.min")
  a=row.names(b_fit_cv)
  b=b_fit_cv[1:length(a)]  
  
  beta_out=data.frame(
    Variables=a[!b==0],
    Beta=b[!b==0],stringsAsFactors = F
  )
  
  list(beta_out=beta_out,lambda=cvfit$lambda.min,prediction=predict(cvfit, newx = x, s = "lambda.min"))
}




######################################## lasso glm, Calculate beta and lambda ########################################

lasso_glm_CalMinLambda=function(indata,var_y,var_x_in,lasso_seed){
  #lasso fit
  indata1=output_nomissing_data(indata,c(var_x_in,var_y))
  N_used=nrow(indata1)
  fit_formula=formula(paste0(" ~ ",paste(var_x_in,collapse = "+")))
  x=model.matrix(fit_formula, indata1)
  y=unlist(indata1[var_y])
  set.seed(lasso_seed)
  cvfit = glmnet::cv.glmnet(x, y,family="gaussian", alpha=1,nfolds = 10)
  
  #extract beta from cv model
  beta_cv=coef(cvfit, cvfit$lambda.min)
  lasso_fit_b1=data.frame(
    Variables=row.names(beta_cv),
    Beta=round(beta_cv[1:length(beta_cv)],4),
    stringsAsFactors = F
  )
  lasso_fit_b1=lasso_fit_b1[!lasso_fit_b1$Beta==0,]
  lasso_fit_b1
  
  #output
  list(
    cvfit=cvfit,
    lambda_out=cvfit$lambda.min,
    Beta_cv=lasso_fit_b1,
    N_used=N_used
  )  
}

######################################## lasso glm, check ########################################


lasso_glm_check=function(indata,var_y,var_x_in,var_adj){
  out_check=bind_rows(lapply(1:1000, function(lasso_seed){
    print(paste0(var_y,": ",lasso_seed))
    lasso_fit_one=lasso_glm_CalMinLambda(indata,var_y,c(var_x_in,var_adj),lasso_seed)
    
    beta_cv=lasso_fit_one$Beta_cv
    beta_cv1=beta_cv[!substr(lasso_fit_one$Beta_cv$Variables,1,6) %in% substr(c("(Intercept)",var_adj),1,6),]
    
    if(nrow(beta_cv1)>0){
      indata$calculated_score=cal_score_basedOnBeta(indata,beta_cv1)
      cor_test=cor.test(indata$calculated_score,unlist(indata[var_y]))
      
      out_one=data.frame(
        random_seed=lasso_seed,
        Outcome=var_y,
        lambda_min=lasso_fit_one$lambda_out,
        N_used=lasso_fit_one$N_used,
        N_beta=nrow(lasso_fit_one$Beta_cv),
        Beta=paste0(paste0(lasso_fit_one$Beta_cv$Variables," (",lasso_fit_one$Beta_cv$Beta,")"),collapse=", "),
        Corr_outcome_calculatedScore=sprintf(paste0("%.",4,"f"),cor_test$estimate),
        P_corr=signif(cor_test$p.value,2),
        stringsAsFactors = F
      )
    }
    
    if(nrow(beta_cv1)==0){
      
      out_one=data.frame(
        random_seed=lasso_seed,
        Outcome=var_y,
        lambda_min=lasso_fit_one$lambda_out,
        N_used=lasso_fit_one$N_used,
        N_beta=nrow(lasso_fit_one$Beta_cv),
        Beta=paste0(paste0(lasso_fit_one$Beta_cv$Variables," (",lasso_fit_one$Beta_cv$Beta,")"),collapse=", "),
        Corr_outcome_calculatedScore=NA,
        P_corr=NA,
        stringsAsFactors = F
      )
    }
    out_one
    
  }))
  out_check
  
}

######################################## lasso logistic, calculate lambda, beta and predicted value ########################################
#Based on
#http://www.sthda.com/english/articles/36-classification-methods-essentials/149-penalized-logistic-regression-essentials-in-r-ridge-lasso-and-elastic-net/

#var_x_in=c(meta_all)
#var_y="traj3_2g"
#var_id="studyid"
#lasso_seed=1000

#names(indata)[1]


#x1=lasso_logistic_CalMinLambda(indata,var_y,var_x_in,"studyid",10)


lasso_logistic_CalMinLambda=function(indata,var_y,var_x_in,var_id,lasso_seed){
  #lasso fit
  indata1=output_nomissing_data(indata,c(var_id,var_x_in,var_y))
  N_used=nrow(indata1)
  fit_formula=formula(paste0(" ~ ",paste(var_x_in,collapse = "+")))
  x=model.matrix(fit_formula, indata1)
  y=as.factor(unlist(indata1[var_y]))
  set.seed(lasso_seed)
  cvfit = glmnet::cv.glmnet(x, y,family="binomial", alpha=1,nfolds = 10)
  
  #extract beta from cv model
  beta_cv=coef(cvfit, cvfit$lambda.min)
  lasso_fit_b1=data.frame(
    Variables=row.names(beta_cv),
    Beta=round(beta_cv[1:length(beta_cv)],4),
    stringsAsFactors = F
  )
  lasso_fit_b1=lasso_fit_b1[!lasso_fit_b1$Beta==0,]
  lasso_fit_b1
  
  #calculate predicted value
  set.seed(lasso_seed)
  fit_lasso_lambdamin=glmnet::glmnet(x, y, alpha = 1, family = "binomial",lambda =cvfit$lambda.min)
  predicted=predict(fit_lasso_lambdamin,newx = x)
  out_predicted=data.frame(
    id=unlist(indata1[var_id]),
    predicted=predicted[1:length(predicted)],
    stringsAsFactors = F
  )
  names(out_predicted)[1]=var_id
  
  #extract beta from no cv model
  Beta_raw=coef(fit_lasso_lambdamin)
  lasso_fit_b=data.frame(
    Variables=row.names(Beta_raw),
    Beta=round(Beta_raw[1:length(Beta_raw)],4),
    stringsAsFactors = F
  )
  lasso_fit_b=lasso_fit_b[!lasso_fit_b$Beta==0,]
  lasso_fit_b
  
  #output
  list(
    cvfit=cvfit,
    lambda_out=cvfit$lambda.min,
    Beta_cv=lasso_fit_b1,
    Beta_NotCV=lasso_fit_b,
    out_predicted=out_predicted,
    N_used=N_used
  )
}


######################################## lasso logistic, check ########################################

lasso_logistic_check=function(indata,var_y,var_x_in,var_adj,var_id){
  
  out_check=bind_rows(lapply(1:1000, function(lasso_seed){
    print(lasso_seed)
    lasso_fit_one=lasso_logistic_CalMinLambda(indata,var_y,c(var_x_in,var_adj),var_id,lasso_seed)
    indata1=indata %>% left_join(lasso_fit_one$out_predicted)
    beta_cv=lasso_fit_one$Beta_cv
    beta_cv1=beta_cv[!substr(lasso_fit_one$Beta_cv$Variables,1,6) %in% substr(c("(Intercept)" ,var_adj),1,6),]
    
    if(nrow(beta_cv1)>0){
      
      
      indata1$calculated_score=cal_score_basedOnBeta(indata,beta_cv1)
      
      out_one=cbind(
        random_seed=lasso_seed,
        Outcome=var_y,
        lambda_min=lasso_fit_one$lambda_out,
        N_beta=nrow(lasso_fit_one$Beta_cv),
        Beta=paste0(paste0(lasso_fit_one$Beta_cv$Variables," (",lasso_fit_one$Beta_cv$Beta,")"),collapse=", "),
        cal_OR_one(indata1,var_y,"calculated_score",4)[c('N_used','OR','CI','P')],
        calculatedScore_AUC=cal_auc(indata1,var_y,'calculated_score',n_digit=4)$auc_ci,
        LassoCVPredictedValue_AUC=cal_auc(indata1,var_y,'predicted',n_digit=4)$auc_ci
      )
      out_one$Beta=as.character(out_one$Beta)
      out_one$calculatedScore_AUC=as.character(out_one$calculatedScore_AUC)
      out_one$LassoCVPredictedValue_AUC=as.character(out_one$LassoCVPredictedValue_AUC)
    }
    
    if(nrow(beta_cv1)==0){
      out_one=data.frame(
        random_seed=lasso_seed,
        Outcome=var_y,
        lambda_min=lasso_fit_one$lambda_out,
        N_beta=nrow(lasso_fit_one$Beta_cv),
        Beta=paste0(paste0(lasso_fit_one$Beta_cv$Variables," (",lasso_fit_one$Beta_cv$Beta,")"),collapse=", "),
        N_used=lasso_fit_one$N_used,OR=NA,CI=NA,P=NA,
        calculatedScore_AUC=NA,
        LassoCVPredictedValue_AUC=NA,stringsAsFactors = F
      )
      
    }
    
    
    out_one
    
  }))
  out_check
}


######################################## lasso logistic, known lambda, calcualte beta and predicted value ########################################
lasso_logistic_KnownLambda=function(indata,var_y,var_x_in,var_id,lambda_in){
  #calculate predicted value
  indata1=output_nomissing_data(indata,c(var_id,var_x_in,var_y))
  fit_formula=formula(paste0(" ~ ",paste(var_x_in,collapse = "+")))
  x=model.matrix(fit_formula, indata1)
  y=as.factor(unlist(indata1[var_y]))
  
  fit_lasso_lambdamin=glmnet::glmnet(x, y, alpha = 1, family = "binomial",lambda =lambda_in)
  predicted=predict(fit_lasso_lambdamin,newx = x)
  out_predicted=data.frame(
    id=unlist(indata1[var_id]),
    predicted=predicted[1:length(predicted)],
    stringsAsFactors = F
  )
  names(out_predicted)[1]=var_id
  
  #extract beta from no cv model
  Beta_raw=coef(fit_lasso_lambdamin)
  lasso_fit_b=data.frame(
    Variables=row.names(Beta_raw),
    Beta=round(Beta_raw[1:length(Beta_raw)],4),
    stringsAsFactors = F
  )
  lasso_fit_b=lasso_fit_b[!lasso_fit_b$Beta==0,]
  lasso_fit_b
  
  #output
  list(
    lasso_fit=fit_lasso_lambdamin,
    Beta_NotCV=lasso_fit_b,
    out_predicted=out_predicted
  )
}
