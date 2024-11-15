########################################  Convert beta and se to OR values ########################################

convert_beta_to_or=function(x,digit_n){
  x1=sprintf(paste0("%.2f"),exp(as.numeric(word(x,1,sep=" "))))
  x2=sprintf(paste0("%.2f"),exp(as.numeric(gsub("[()]","",word(x,2,sep=" ")))))
  x3=sprintf(paste0("%.2f"),exp(as.numeric(gsub("[()]","",word(x,4,sep=" ")))))

  paste0(x1," (",x2," - ",x3,")")
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

########################################  Test missing ########################################
output_nomissing_data=function(indata,var_all){
  indata1=indata[var_all]
  indata1$missing_index=apply(apply(indata1,2,function(x){ifelse(is.na(x),1,0)}),1,sum)
  outdata=indata[indata1$missing_index==0,]
}

######################################## logistic, one outcome (2 levels), one x or multiple x (with adjustments)  ##################################
#digit_n=2


#outcome_one="traj3_2g"
#var_x=c("score_lasso_traj_bySD",model1)

cal_OR_one=function(indata,outcome_one,var_x,digit_n){
  indata$outcome=factor(unlist(indata[outcome_one]))
  fit_formula=formula(paste0("outcome"," ~ ",paste0(var_x,collapse = "+")))
  logit=glm(fit_formula,family = "binomial",data=indata)

  style_number=paste0("%.",digit_n,"f")

  #SE(logit)

  out_one=data.frame(
    Outcome=outcome_one,
    Variable=var_x[1],
    Cov=paste0(var_x,collapse = "+"),
    N_used=logit$df.null+1,
    var_x_level=names(coef(logit)),
    Beta=sprintf(style_number,coef(logit)),
    OR=sprintf(style_number,exp(coef(logit))),
    CI=paste0(sprintf(style_number,exp(confint.default(logit)[,1])),"-",
              sprintf(style_number,exp(confint.default(logit)[,2]))),
    P_raw=signif(summary(logit)$coefficients[,4],digits = 2),
    P=convert_p(summary(logit)$coefficients[,4]),
    stringsAsFactors = F
  )

  out_one=out_one[-1,]
  out_one
}



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


cal_rr_possion=function(indata,var_outcome_list,var_x_list,adj,digit_n){
  style_number=paste0("%.",digit_n,"f")
  out_logit=bind_rows(lapply(var_outcome_list, function(var_outcome_one){

    bind_rows(lapply(var_x_list, function(var_x_one){
      indata$outcome=as.numeric(factor(unlist(indata[var_outcome_one])))-1

      if(is.na(adj)[1]|adj[1]==""){fit_formula=formula(paste0("outcome"," ~ ",paste0(c(var_x_one),collapse = "+")))}
      if(!(is.na(adj)[1]|adj[1]=="")){fit_formula=formula(paste0("outcome"," ~ ",paste0(c(var_x_one,adj),collapse = "+")))}

      logit=glm(fit_formula, family="poisson", data=indata)

      out_one=data.frame(
        Outcome=var_outcome_one,
        N_used=logit$df.null+1,
        Cov=paste0(adj,collapse = "+"),
        Var_x=var_x_one,
        x_level=substr(names(coef(logit)),nchar(var_x_one)+1,nchar(names(coef(logit)))),
        Beta=sprintf(style_number,coef(logit)),
        RR=sprintf(style_number,exp(coef(logit))),
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


######################################## logistic, one or multiple outcome, one categorical x, with trend p values  ##################################


cal_or_categorical_with_trend=function(indata,var_y_list,var_x_one_cat,var_x_one_con,digit_n){
  bind_rows(lapply(var_y_list, function(var_y_one){
    temp1=cal_freq_NoGroup(indata,var_x_one_cat,digit_n)
    temp2=cal_or(indata,var_y_one,var_x_one_cat,NA,digit_n)
    temp3=cal_or(indata,var_y_one,var_x_one_con,NA,digit_n)

    temp2$N_used=NULL

    names(temp1)=c('Var_x',"x_level","N_used","Level_N")

    data_out=temp1 %>% left_join(temp2)

    data_out$Beta=0
    data_out$OR[2]=1
    data_out$CI[2]="Reference"

    data_out$OR_CI=paste0(data_out$OR," (",data_out$CI,")")
    data_out$Outcome=var_y_one
    data_out$N_used=data_out$N_used[2]

    data_out=data_out[2:nrow(data_out),]

    data_out1=dcast(
      data_out,
      Outcome+N_used~x_level,
      value.var="OR_CI")

    data_out1$P_raw=temp3$P_raw

    data_out1$P_trend=temp3$P
    data_out1
  }))
}


######################################## genenate predicted value  ##################################

#indata
#outcome
#var_in

Cal_predicted=function(indata,outcome,var_in){
  indata$outcome=factor(unlist(indata[outcome]))
  indata$outcome=relevel(indata$outcome,ref =levels(indata$outcome)[1])
  formula_fit=formula(paste0('outcome',"~",paste0(var_in,collapse = "+")))
  fit_logit=glm(formula_fit,family = "binomial",data=indata)
  predict(fit_logit)
}
######################################## calculate auc  ##################################
#cal_auc(data_train,"OS_status","score")

#indata=data_for_predict
#outcome="obesity_2level"
#var_in=cov_model3


cal_auc=function(indata,outcome,var_in,n_digit=2){
  indata=output_nomissing_data(indata,c(outcome,var_in))

  indata$outcome=factor(unlist(indata[outcome]))
  indata$outcome=relevel(indata$outcome,ref =levels(indata$outcome)[1])

  formula_fit=formula(paste0('outcome',"~",paste0(var_in,collapse = "+")))
  fit_logit=glm(formula_fit,family = "binomial",data=indata)

  fit_roc=pROC::roc(unlist(indata[outcome]),predict(fit_logit),ci=T)

  auc_ci=paste0(sprintf(paste0("%.",n_digit,"f"),fit_roc$auc)," (",
                paste0(sprintf(paste0("%.",n_digit,"f"),fit_roc$ci[c(1,3)]),collapse = " - "),")")

  youdeng_raw=fit_roc$sensitivities+fit_roc$specificities-1

  youdeng=max(youdeng_raw)

  sensitivity=round(fit_roc$sensitivities[youdeng_raw==youdeng],4)
  specificity=round(fit_roc$specificities[youdeng_raw==youdeng],4)

  data.frame(
    Outcome=outcome,
    Model=paste0(var_in,collapse = "+"),
    auc_ci=auc_ci,
    stringsAsFactors = F

    #,
    #youden=youdeng,
    #sensitivity=sensitivity,
    #specificity=specificity
  )
}

########################################### Calculate score based on beta values ##################################################

#data_beta should have tow columns, the name of variables(Variables) and the beta values(Beta).

cal_score_basedOnBeta=function(indata,data_beta){
  score=0
  for (j in 1:length(data_beta$Variables)){
    score=score+unlist(indata[data_beta$Variables[j]])*(data_beta$Beta[j])}
  score
}


######################################## calculate score based on beta values ########################################

cal_score_new=function(in_data,in_beta){
  in_data1=in_data %>% mutate(hp=as.numeric(substr(hp,1,1)),
                              DM=as.numeric(substr(DM,1,1)),
                              dyslipidemia=as.numeric(substr(dyslipidemia,1,1)),
                              physical_activity=as.numeric(substr(physical_activity,1,1)),
                              favor_salt=as.numeric(substr(favor_salt,1,1)),
                              chd_family=as.numeric(substr(chd_family,1,1)),
                              stroke_family=as.numeric(substr(stroke_family,1,1))

  )

  in_data1$physical_activity=ifelse(in_data1$physical_activity==2,1,0)

  print(bind_rows(lapply(c("hp","DM","dyslipidemia","physical_activity","favor_salt","chd_family","stroke_family"), function(var_one){
    cal_freq_NoGroup(in_data1,var_one,2)
  }))
  )

  in_beta1=in_beta %>%
    filter(!Variable=="Intercept") %>%
    filter(!(Variable=="physical_activity" & substr(word(level,2,sep=" "),1,1)==1)) %>%
    mutate(beta=as.numeric(word(beta_se,1,sep=" ")))

  for(i in 1:nrow(in_beta1)){
    print(i)
    value_one=(in_beta1$beta[i]*unlist(in_data1[in_beta1$Variable[i]]))
    if(i==1){value=value_one}
    if(i> 1){value=value+value_one}
  }

  data.frame(
    obs=in_data1$obs,
    ScoreCVDNew=round(value,4),stringsAsFactors = F
  )
}

######################################## calculate HL test ########################################


#install.packages("ResourceSelection")

#ResourceSelection::hoslem.test(score_male$cvd1, score_male$ScoreCVDNew, g=10)
#ResourceSelection::hoslem.test(score_male$cvd1, score_male$Score_cvd, g=10)

hosmerlem = function(y, yhat, g=10) {
  cutyhat = cut(yhat,
                breaks = quantile(yhat, probs=seq(0,
                                                  1, 1/g)), include.lowest=TRUE)
  obs = xtabs(cbind(1 - y, y) ~ cutyhat)
  expect = xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq = sum((obs - expect)^2/expect)
  P = 1 - pchisq(chisq, g - 2)
  return(list(chisq=chisq,p.value=P))
}

#logreg = glm(homeless ~ female + i1 + cesd + age + substance,family=binomial)
#hosmerlem(y=homeless, yhat=fitted(logreg))


######################################## Auc compare ########################################
#indata=data_for_predict
#outcome="obesity_2level"
#var_base=cov_model3
#var_compare=c(cov_model3,"score_obesity")

#indata=data_for_auc
#outcome="CV4Obesity_2g"
#var_base=model1
#var_compare=c(model1,"score_lasso_obesity")
#digit_n=3

auc_compare=function(indata,outcome,var_base,var_compare,digit_n=1){
  var_in=unique(c(outcome,var_base,var_compare))
  indata=output_nomissing_data(indata,var_in)

  indata$outcome=factor(unlist(indata[outcome]))
  indata$outcome=relevel(indata$outcome,ref =levels(indata$outcome)[1])

  formula_fit_base=   formula(paste0('outcome',"~",paste0(var_base,collapse = "+")))
  formula_fit_compare=formula(paste0('outcome',"~",paste0(var_compare,collapse = "+")))

  fit_logit_base=glm(formula_fit_base,family = "binomial",data=indata)
  fit_logit_compare=glm(formula_fit_compare,family = "binomial",data=indata)

  fit_roc_base=pROC::roc(unlist(indata[outcome]),predict(fit_logit_base),ci=T)
  fit_roc_compare=pROC::roc(unlist(indata[outcome]),predict(fit_logit_compare),ci=T)

  out_auc_compare=pROC::roc.test(fit_roc_base,fit_roc_compare)

  #calculate auc difference.
  AUC_dif=out_auc_compare$roc2$auc-out_auc_compare$roc1$auc
  AUC_dif_se=abs(AUC_dif/out_auc_compare$statistic)
  AUC_dif_ul=sprintf(paste0("%.",digit_n,"f"),AUC_dif+qnorm(1-0.05/2)*AUC_dif_se)
  AUC_dif_ll=sprintf(paste0("%.",digit_n,"f"),AUC_dif-qnorm(1-0.05/2)*AUC_dif_se)

  #output
  data.frame(
    outcome=outcome,
    var_base=paste0(var_base,collapse = "+"),
    var_compare=paste0(var_compare,collapse = "+"),
    AUC_base=paste0(
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc1$auc)," (",
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc1$ci[1]),"-",
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc1$ci[3]),")"),

    AUC_for_compare=paste0(
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc2$auc)," (",
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc2$ci[1]),"-",
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc2$ci[3]),")"),

    AUC_difference_ci=paste0(sprintf(paste0("%.",digit_n,"f"),AUC_dif)," (",
                      AUC_dif_ll,"-",AUC_dif_ul,")"),
    Z_statistics=out_auc_compare$statistic,
    P_auc_compare_raw=signif(out_auc_compare$p.value,2),
    P_auc_compare=convert_p(out_auc_compare$p.value),
    stringsAsFactors = F
  )
}



auc_compare_usingVector=function(outcome_base,var_base,outcome_compare,var_compare,digit_n=3){
  df_base=data.frame(  outcome=as.factor(outcome_base),  var=var_base,  stringsAsFactors = F) %>% drop_na()
  df_compare=data.frame(  outcome=as.factor(outcome_compare),  var=var_compare,  stringsAsFactors = F) %>% drop_na()

  df_base$outcome=relevel(df_base$outcome,ref =levels(df_base$outcome)[1])
  df_compare$outcome=relevel(df_compare$outcome,ref =levels(df_compare$outcome)[1])

  formula_fit_base=   formula(paste0('outcome',"~",paste0(var_base,collapse = "+")))
  formula_fit_compare=formula(paste0('outcome',"~",paste0(var_compare,collapse = "+")))

  fit_logit_base=glm(outcome~var,family = "binomial",data=df_base)
  fit_logit_compare=glm(outcome~var,family = "binomial",data=df_compare)

  fit_roc_base=pROC::roc(unlist(df_base$outcome),predict(fit_logit_base),ci=T)
  fit_roc_compare=pROC::roc(unlist(df_compare$outcome),predict(fit_logit_compare),ci=T)

  out_auc_compare=pROC::roc.test(fit_roc_base,fit_roc_compare)

  #calculate auc difference.
  AUC_dif=out_auc_compare$roc2$auc-out_auc_compare$roc1$auc
  AUC_dif_se=abs(AUC_dif/out_auc_compare$statistic)
  AUC_dif_ul=sprintf(paste0("%.",digit_n,"f"),AUC_dif+qnorm(1-0.05/2)*AUC_dif_se)
  AUC_dif_ll=sprintf(paste0("%.",digit_n,"f"),AUC_dif-qnorm(1-0.05/2)*AUC_dif_se)

  #output
  data.frame(
    AUC_base=paste0(
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc1$auc)," (",
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc1$ci[1]),"-",
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc1$ci[3]),")"),

    AUC_for_compare=paste0(
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc2$auc)," (",
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc2$ci[1]),"-",
      sprintf(paste0("%.",digit_n,"f"),out_auc_compare$roc2$ci[3]),")"),

    AUC_difference_ci=paste0(sprintf(paste0("%.",digit_n,"f"),AUC_dif)," (",
                             AUC_dif_ll,"-",AUC_dif_ul,")"),
    Z_statistics=out_auc_compare$statistic,
    P_auc_compare_raw=signif(out_auc_compare$p.value,2),
    P_auc_compare=convert_p(out_auc_compare$p.value),
    stringsAsFactors = F
  )
}


#get_auc_bootstrap ------------------------
get_auc_bootstrap=function(indata,var_y,var_x,times=1000){
  indata$y=unlist(indata[var_y])
  indata$x=unlist(indata[var_x])

  auc_boot=rep(NA,times)
  index=1:nrow(indata)
  for(i in 1:times){
    # print(i)
    set.seed(i)
    index_i=sample(index,replace = T)
    auc_boot[i]=auc(indata$y[index_i],indata$x[index_i])
  }

  data.frame(var_y=var_y,var_x=var_x,auc=auc_boot,stringsAsFactors = F)
}

# df_in=df_auc_bootstrap
# var_x1="PA_score"
# var_x2="Tumor_differentiation"

compare_auc=function(df_in,var_x1,var_x2){
  value1=unlist(df_in$auc[df_in$var_x==var_x1],)
  value2=unlist(df_in$auc[df_in$var_x==var_x2],)

  t_out=t.test(value1,value2)

  df_out=data.frame(var_x1=var_x1,var_x2=var_x2,
                    t=t_out$statistic,P_raw=t_out$p.value,
                    stringsAsFactors = F)
  df_out
}


######################################## Auc plot ########################################

#auc_plot(data_male_with_score,"cvd","age","Age","Male",NA)

#auc_plot(data_male_with_score,"cvd",list("score","age"),c("Score","Age"),"Male",2)


auc_plot=function(indata,outcome,var_list,label_list,plot_label,ref){
  color_in=RColorBrewer::brewer.pal(n = 9, name = 'Set1')[1:length(var_list)]

  indata=output_nomissing_data(indata,c(outcome,unlist(var_list)))

  indata$outcome=factor(unlist(indata[outcome]))
  indata$outcome=relevel(indata$outcome,ref =levels(indata$outcome)[1])

  cal_plot_para=function(var_in){
    print(paste0("var = ",paste0(var_in,collapse = "+")))
    formula_fit=formula(paste0('outcome',"~",paste0(var_in,collapse = "+")))

    fit_logit=glm(formula_fit,family = "binomial",data=indata)

    fit_roc=pROC::roc(unlist(indata[outcome]),predict(fit_logit),ci=T)

    out=data.frame(
      spec=fit_roc$specificities,
      sens=fit_roc$sensitivities,stringsAsFactors = F
    )
    out
  }

  plot(NA,type="l",xlim=c(0,1),ylim=c(0,1),
       xlab = "1 - Specificity",ylab = "Sensitivity",main=plot_label)

  abline(h = seq(0,1,0.1),col = "lightgray",lty = 3)
  abline(v = seq(0,1,0.1),col = "lightgray",lty = 3)

  if(!is.na(ref)){
    plot_para_one=cal_plot_para(var_list[[ref]])
    polygon(x = 1-plot_para_one$spec,y = plot_para_one$sens,col = "lightgray")
    polygon(x = c(0,1,1),y = c(0,0,1),col = "lightgray")
  }

  abline(a=0,b=1,col="gray65")
  for(i in 1:length(var_list)){
    plot_para_one=cal_plot_para(var_list[[i]])
    lines(1-plot_para_one$spec,plot_para_one$sens,lwd=2,col=color_in[i])
  }
  legend("bottomright",
         col=color_in,
         lwd=2,bty="n",
         legend=label_list)
}




######################################## Auc plot for Qi ########################################
#indata=data_train
#outcome="Recurrence_status"
#var_list=list("predicted","stage_1","stage_2","stage_3",'stage_4')
#label_list=c("Predicted","S1","S2","S3",'S4')
#var_in=var_list[2]

#set.seed(5)
#color_in=sample(colors(),length(var_list))


auc_plot_Qi=function(indata,outcome,var_list,label_list,plot_label,ref){
  indata=output_nomissing_data(indata,c(outcome,unlist(var_list)))

  indata$outcome=factor(unlist(indata[outcome]))
  indata$outcome=relevel(indata$outcome,ref =levels(indata$outcome)[1])

  cal_plot_para=function(var_in){
    print(paste0("var = ",paste0(var_in,collapse = "+")))
    formula_fit=formula(paste0('outcome',"~",paste0(var_in,collapse = "+")))

    fit_logit=glm(formula_fit,family = "binomial",data=indata)

    fit_roc=pROC::roc(unlist(indata[outcome]),predict(fit_logit),ci=T)

    out=data.frame(
      spec=fit_roc$specificities,
      sens=fit_roc$sensitivities,stringsAsFactors = F
    )
    out
  }

  plot(NA,type="l",xlim=c(0,1),ylim=c(0,1),
       xlab = "1 - Specificity",ylab = "Sensitivity",main=plot_label)

  abline(h = seq(0,1,0.1),col = "lightgray",lty = 3)
  abline(v = seq(0,1,0.1),col = "lightgray",lty = 3)

  if(!is.na(ref)){
    plot_para_one=cal_plot_para(var_list[[ref]])
    polygon(x = 1-plot_para_one$spec,y = plot_para_one$sens,col = "lightgray")
    polygon(x = c(0,1,1),y = c(0,0,1),col = "lightgray")
  }

  abline(a=0,b=1,col="gray65")
  for(i in 1:length(var_list)){
    plot_para_one=cal_plot_para(var_list[[i]])
    lines(1-plot_para_one$spec,plot_para_one$sens,lwd=2,col=color_in[i])
  }
  legend(x = 0.5,y=0.2,
         col=color_in,
         lwd=2,bty="n",
         legend=label_list)
}



