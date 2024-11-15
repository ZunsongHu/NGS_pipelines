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


########################################### Calculate score based on beta values ##################################################
cal_score_basedOnBeta=function(indata,data_beta){
  score=0
  for (j in 1:length(data_beta$Variables)){
    score=score+unlist(indata[data_beta$Variables[j]])*(data_beta$Beta[j])}
  score
}


########################################### cox output all variable results ##################################################
cal_cox_output_multiple=function(indata,outcome,outcome_t,var_in,digit_n){

  fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var_in,collapse = "+")))
  cox_fit=survival::coxph(fit_formula,data=indata)
  style_number=paste0("%.",digit_n,"f")

  out_one=data.frame(
    Outcome=outcome,
    Adj=ifelse(is.na(var_in)[1],NA,paste0(c(var_in),collapse = "+")),
    N_used=cox_fit$n,

    Variable_level=names(cox_fit$coefficients),

    Beta_se=paste0(sprintf(style_number,cox_fit$coefficients)," (",
                   sprintf(style_number,(cox_fit$coefficients-confint(cox_fit)[1])/1.96),
                   ")"),
    RR_CI=paste0(sprintf(style_number,exp(cox_fit$coefficients))," (",
                 sprintf(style_number,exp(confint(cox_fit)[,1])),"-",
                 sprintf(style_number,exp(confint(cox_fit)[,2])),")"),

    P_raw=signif(summary(cox_fit)$coefficients[,5],digits = 3),
    P=convert_p(summary(cox_fit)$coefficients[,5]),
    stringsAsFactors = F
  )
  out_one
}


########################################### cox, one or multiple x, one or multiple y ##################################################
#indata=data_pheno_withLassoScore1[data_pheno_withLassoScore1$datasource=='test1',]
#outcome="pfs.status"
#outcome_t="pfs"
#var_x="lasso_score_summed"
#digit_n=3
#var_adj=adj_0.25

library(survival)
cal_cox=function(indata,outcome_list,outcome_t,var_x_list,var_adj,digit_n){
  style_number=paste0("%.",digit_n,"f")

  bind_rows(lapply(outcome_list, function(outcome){
    bind_rows(lapply(var_x_list, function(var_x){

      if(is.na(var_adj)[1]){indata1=indata[c(outcome,outcome_t,var_x)] %>% na.omit}
      if(!is.na(var_adj)[1]){indata1=indata[c(outcome,outcome_t,var_x,var_adj)] %>% na.omit}

      #print(nrow(indata1))
      if(nrow(indata1)>1 & !any(apply(indata1,2,function(x){length(levels(as.factor(x)))})==1)){
        indata1[,1]=as.numeric(as.factor(unlist(indata1[1])))

        var_x_in=ifelse(is.na(var_adj)[1],var_x,paste0(c(var_x,var_adj),collapse = "+"))

        fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",var_x_in))
        cox_fit=survival::coxph(fit_formula,data=indata1)

        #summarytools::view(summarytools::dfSummary(indata1))

        out_one=data.frame(
          Outcome=outcome,
          Variable=substr(names(cox_fit$coefficients),1,nchar(var_x)),
          Adj=ifelse(is.na(var_adj)[1],NA,paste0(c(var_adj),collapse = "+")),
          N_used=cox_fit$n,

          Var_level=substr(names(cox_fit$coefficients),(nchar(var_x)+1),nchar(names(cox_fit$coefficients))),

          Beta_se=paste0(sprintf(style_number,cox_fit$coefficients)," (",
                         sprintf(style_number,(cox_fit$coefficients-confint(cox_fit)[1])/1.96),
                         ")"),
          RR_CI=paste0(sprintf(style_number,exp(cox_fit$coefficients))," (",
                       sprintf(style_number,exp(confint(cox_fit)))[1],"-",
                       sprintf(style_number,exp(confint(cox_fit)))[2],")"),

          P_raw=signif(summary(cox_fit)$coefficients[,5],digits = 3),
          P=convert_p(summary(cox_fit)$coefficients[,5]),
          stringsAsFactors = F
        )
        out_one=out_one[out_one$Variable==var_x,]
      }
      if(nrow(indata1)<=1 | any(apply(indata1,2,function(x){length(levels(as.factor(x)))})==1)){
        out_one=data.frame(
          Outcome=outcome,
          Variable=paste0(var_x,collapse = "+"),
          Adj=ifelse(is.na(var_adj)[1],NA,paste0(c(var_adj),collapse = "+")),
          N_used=0,
          Var_level=NA,
          Beta_se=NA,
          RR_CI=NA,
          P_raw=NA,
          P=NA,
          stringsAsFactors = F
        )
      }
      out_one

    }))
  }))
}


########################################### cox_multiple_var_customed ##################################################
#indata=data_test_new
#outcome="Recurrence_status"
#outcome_t="RFS"
#round_number=4
#var="V4_g"

#var_all=var_all_train

cox_multiple_var_customed=function(indata,outcome,outcome_t,var_all){
  indata[var_all]=apply(indata[var_all], 2,factor)

  out_levels=bind_rows(lapply(var_all,function(var){
    freq_one=data.frame(table(indata[var]),stringsAsFactors = F)
    names(freq_one)=c("Predictor","Freq")
    freq_one=as.data.frame(apply(freq_one,2, as.character),stringsAsFactors = F)
    freq_one$label=paste0(var,freq_one$Predictor)
    rbind(
      data.frame(Predictor=var,Freq=NA,label=NA,stringsAsFactors = F),
      freq_one)
  }))


  fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var_all,collapse = "+")))
  cox_fit=coxph(fit_formula, data = indata)

  cox_fit_out=data.frame(
    label=names(cox_fit$coefficients),
    rr=paste0(round(exp(cox_fit$coefficients),round_number)," (",
              round(exp(confint(cox_fit)),round_number)[,1],"-",
              round(exp(confint(cox_fit)),round_number)[,2],")"),
    p=signif(summary(cox_fit)$coefficients[,5],digits = 2),
    Beta=round(coef(cox_fit),round_number),
    stringsAsFactors = F)

  rr_out=out_levels %>% left_join(cox_fit_out)
  rr_out$rr=ifelse(rr_out$Predictor=="0" & is.na(rr_out$rr),'1 [Reference]',rr_out$rr)
  rr_out$Points=ifelse(rr_out$Predictor %in% c("0",'1','2',"3"),rr_out$Predictor,NA)
  rr_out
}


######################################## cox, one or multiple outcome, one categorical x, with trend p values  ##################################


cal_cox_categorical_with_trend=function(indata,outcome_list,outcome_t,var_x,var_adj,digit_n){
  bind_rows(
    lapply(outcome_list, function(outcome){
      style_number=paste0("%.",digit_n,"f")

      if(is.na(var_adj)){indata1=indata[c(outcome,outcome_t,var_x)] %>% na.omit}
      if(!is.na(var_adj)){indata1=indata[c(outcome,outcome_t,var_x,var_adj)] %>% na.omit}

      indata1[,1]=as.numeric(as.factor(indata1[,1]))

      indata1$var_x=as.numeric(as.factor(unlist(indata1[var_x])))


      freq_one=data.frame(table(indata1[var_x]),stringsAsFactors = F)
      names(freq_one)=c("Var_level","Freq")
      freq_one=as.data.frame(apply(freq_one,2, as.character),stringsAsFactors = F)
      freq_one$label=paste0(var_x,"_",freq_one$Var_level," (N=",freq_one$Freq,")")


      var_x_in=ifelse(is.na(var_adj),var_x,paste0(c(var_x,var_adj),collapse = "+"))
      fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",var_x_in))
      cox_fit=coxph(fit_formula,data=indata1)


      var_x_in1=ifelse(is.na(var_adj),"var_x",paste0(c("var_x",var_adj),collapse = "+"))
      fit_formula1=formula(paste0("Surv(",outcome_t,",",outcome,")~",var_x_in1))
      cox_fit1=coxph(fit_formula1,data=indata1)


      out_one=data.frame(
        Variable=substr(names(cox_fit$coefficients),1,nchar(var_x)),

        Var_level=substr(names(cox_fit$coefficients),(nchar(var_x)+1),nchar(names(cox_fit$coefficients))),

        RR_CI=paste0(sprintf(style_number,exp(cox_fit$coefficients))," (",
                     sprintf(style_number,exp(confint(cox_fit)))[1],"-",
                     sprintf(style_number,exp(confint(cox_fit)))[2],")"),

        P_raw=signif(summary(cox_fit)$coefficients[,5],digits = 3),
        P=convert_p(summary(cox_fit)$coefficients[,5]),
        stringsAsFactors = F
      )
      out_rr=out_one[out_one$Variable==var_x,]

      out_rr=freq_one %>% left_join(out_rr)
      out_rr$RR_CI[1]="1 (Ref)"
      out_rr$N_used=cox_fit$n
      out_rr$Outcome=outcome
      out_rr$Variable=var_x
      out_rr$Adj=ifelse(is.na(var_adj),NA,paste0(c(var_adj),collapse = "+"))

      out_rr1=dcast(
        out_rr, #name of dataset
        Outcome+Variable+Adj+N_used~label,
        function(x){paste0(unique(x),collapse = "+")},
        value.var="RR_CI")

      out_one1=data.frame(
        Variable=substr(names(cox_fit1$coefficients),1,nchar(var_x)),


        P_raw=signif(summary(cox_fit1)$coefficients[,5],digits = 3),
        P=convert_p(summary(cox_fit1)$coefficients[,5]),
        stringsAsFactors = F
      )
      out_one1=out_one1[out_one1$Variable=='var_x',]
      out_rr1$P_trend_raw=out_one1$P_raw
      out_rr1$P_trend=out_one1$P
      out_rr1
    }))
}


########################################### cox_stepwise ##################################################


cox_stepwise=function(indata,outcome,outcome_t,var_x_all,stepMethod){

style_number=paste0("%.",digit_n,"f")


fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var_x_all,collapse = "+")))
cox_fit=coxph(fit_formula, data = indata)

#sink(paste0(outdir,"stepwise.txt"))
cox_fit_step=step(cox_fit,scope = list(upper=cox_fit),direction = "forward")

cox_fit=cox_fit_step

out_one=data.frame(
  Outcome=outcome,
  Variables=names(coef(cox_fit)),
  Beta_se=paste0(sprintf(style_number,cox_fit$coefficients)," (",
                 sprintf(style_number,(cox_fit$coefficients-confint(cox_fit)[1])/1.96),
                 ")")   ,
  RR_CI=paste0(sprintf(style_number,exp(cox_fit$coefficients))," (",
               sprintf(style_number,exp(confint(cox_fit)))[1],"-",
               sprintf(style_number,exp(confint(cox_fit)))[2],")")
  ,

  P_raw=signif(summary(cox_fit)$coefficients[,5],digits = 3),
  P=convert_p(summary(cox_fit)$coefficients[,5]),
  stringsAsFactors = F
)



out_one
}


#cox_fit=stepwise.coxph_new(Time = outcome_t, Status = outcome, variable.list = var_x_all,sle = 0.1, sls = 0.1,data = indata)

########################################### Calculate predicted value ##################################################
#indata=data_train
#outcome="Recurrence_status"
#outcome_t="RFS"
#var_in=c("V4",'V7','V10','V14','V15','V17','V18')

cal_predicted=function(indata,outcome,outcome_t,var_in){
  fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var_in,collapse = "+")))
  cox_fit=coxph(fit_formula, data = indata)
  predict(cox_fit)
}

########################################### Calculate predicted risk at specified times ##################################################
library(rms)
library(riskRegression)

#indata=data_train_withScoreWeighted
#var="points_sum"
#cal_predicted(data_train_new,'Recurrence_status','RFS',"points_sum")

cal_predicted_year=function(indata,outcome,outcome_t,var){
  ddist=datadist(indata)
  options(datadist='ddist')
  fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var,collapse = "+")))
  cox_fit=cph(fit_formula,data=indata,x=TRUE, y=TRUE,surv = T,time.inc =12)
  risk_predict=data.frame(riskRegression::predictRisk(cox_fit,newdata=indata,times=c(12,36,60)))
  names(risk_predict)=c("predicted_1Year","predicted_3Year","predicted_5Year")
  cbind(indata,risk_predict)
}


########################################### c_index ##################################################


cal_c=function(indata,outcome,outcome_t,var_c,fit_name,digit_n){
  indata1=na.omit(indata[c(outcome,outcome_t,var_c)])
  if(nrow(indata1)<=1){
    c_result=data.frame(fit_name=fit_name,
                        c=NA,
                        ll=NA,
                        ul=NA,
                        c_ci=NA,
                        stringsAsFactors = F)
  }
  if(nrow(indata1)>1){
    c_result=as.data.frame(matrix(NA,nrow=1,ncol=5))
    names(c_result)=c('fit_name','c','ll','ul','c_ci')
    c_result$fit_name=fit_name
    extract_c=function(x){
      c_index=summary(x)
      c_result=data.frame(fit_name=fit_name,
                          c=round(c_index$concordance[1],digit_n),
                          ll=round(c_index$concordance[1]-1.96*c_index$concordance[2],digit_n),
                          ul=round(c_index$concordance[1]+1.96*c_index$concordance[2],digit_n),
                          c_ci=paste0(signif(c_index$concordance[1],digits = 2)," (",
                                      round((c_index$concordance[1]-1.96*c_index$concordance[2]),digit_n)," - ",
                                      round((c_index$concordance[1]+1.96*c_index$concordance[2]),digit_n),")"),
                          stringsAsFactors = F)

      c_result
    }

    if(length(var_c)==1){
      if(nrow(as.data.frame(table(unlist(indata[var_c]))))>1){
        fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var_c,collapse = "+")))
        cox_fit=coxph(fit_formula, data = indata)
        c_result=extract_c(cox_fit)
      }
    }

    if(length(var_c)>1){
      fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var_c,collapse = "+")))
      cox_fit=coxph(fit_formula, data = indata)
      c_result=extract_c(cox_fit)
    }
  }
  c_result
}


########################################### C index compare ##################################################

c_compare=function(indata,outcome,outcome_t,var_base,var_compare,label){
  fit_formula_base=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var_base,collapse = "+")))
  fit_formula_compare=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var_compare,collapse = "+")))

  predict_base=predict(coxph(fit_formula_base, data = indata))
  predict_compare=predict(coxph(fit_formula_compare, data = indata))

  data.frame(
    label=label,
    Model_base=paste0(var_base,collapse = "+"),
    Model_compare=paste0(var_compare,collapse = "+"),
    p_compare=signif(compareC::compareC(unlist(indata[outcome_t]), unlist(indata[outcome]),
                                        unlist(predict_base),unlist(predict_compare))$pval,
                     digits = 2),
    stringsAsFactors = F)

}


########################################### c index bootstrap ##################################################

#indata=data_predict
#outcome="OS_censor"
#outcome_t="OS"
#var="Model2"
#group=4
#bootstrap_time=50

#out_bootstrap_model1=Cal_c_HL_bootstrap_o(data_predict,"OS_censor","OS","Model1",4,1000)

library(survival)
library(rms)

Cal_c_HL_bootstrap_o=function(indata,outcome,outcome_t,var,group,bootstrap_time){
  cal_c=function(x){
    fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var,collapse = "+")))
    summary(coxph(fit_formula, data = x))$concordance[1]
  }

  cal_HL=function(x){
    ddist=datadist(x)
    options(datadist='ddist')
    fit_formula=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var,collapse = "+")))
    cox_multi3=cph(fit_formula,data=x,x=TRUE, y=TRUE,surv = T,time.inc =12*5)
    n_g=floor(nrow(x)/group)
    sink("a.txt")
    cal=calibrate(cox_multi3, cmethod = 'KM',method = 'boot', u = 12*5, m = n_g, B = 1)
    sink()
    pred=cal[,7];obse=cal[,8];se=cal[,10]
    HL_test=data.frame(t(
      data.frame(c(pred,obse,se,MKmisc::HLgof.test(pred3,obse3,group)$C$statistic,MKmisc::HLgof.test(pred3,obse3,group)$C$p.value))
    ))
    row.names(HL_test)=NULL
    names(HL_test)=c(
      'pred1','pred2','pred3','pred4',
      'obse1','obse2','obse3','obse4',
      'se1','se2','se3','se4',
      'statistic_HL','p_value_HL'
    )
    return(HL_test)
  }

  Cal_c_HL_bootstrap=function(i){
    set.seed(i)
    id_sampled=sample(indata$id,replace = T)
    data_one=indata[match(id_sampled,indata$id),]
    out=cal_HL(data_one)
    out$c=cal_c(data_one)
    out
  }

  out_bootstrap=bind_rows(lapply(1:bootstrap_time, Cal_c_HL_bootstrap),.id = "Seed_number")
  list(out_bootstrap=out_bootstrap,
       c_boot=data.frame(c_boot=round(mean(out_bootstrap$c),2)))
}


########################################### calculate cutoff ##################################################
#indata=data_train
#outcome='Recurrence_status'
#outcome_t='RFS'
#var_in="V7"

#data_train["V7"]

#cal_cutoff_cox(data_train,"OS_status","OS",var_in)

cal_cutoff_cox=function(indata,outcome,outcome_t,var_in){
  formula_fit=formula(paste0("Surv(",outcome_t,",",outcome,")~",var_in))
  cox_fit=coxph(formula_fit,data = indata)
  cutoff=survMisc::cutp(cox_fit)
  value_cutoff=unlist(cutoff[[1]][1,1])
  out=data.frame(var_name=var_in,
                 cutoff=value_cutoff,stringsAsFactors = F)
  out
}



########################################### calculate NRI ##################################################
#indata=data_test_with_predicted_Risk
#predicted_value="predicted_1Year"
#outcome="Recurrence_status"
#outcome_t="RFS"

#var_std="stage_4"
#var_new=group

cal_nri=function(indata,outcome,outcome_t,var_std,var_new,predicted_value,label){
  formula_std=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var_std,collapse = "+")))
  formula_new=formula(paste0("Surv(",outcome_t,",",outcome,")~",paste0(var_new,collapse = "+")))

  mstd = coxph(formula_std, data=indata, x=TRUE)
  mnew = coxph(formula_new, data=indata, x=TRUE)

  low=quantile(unlist(indata[predicted_value]),0.1)
  high=quantile(unlist(indata[predicted_value]),0.6)

  aa=nricens::nricens(mdl.std = mstd, mdl.new = mnew, t0 = 12, cut = c(low,high),niter = 10)

  data.frame(
    label=label,
    var_std=var_std,
    var_new=var_new,
    NRI=aa$nri[1,1],
    P_NRI=aa$nri[4,1],stringsAsFactors = F
  )
}
