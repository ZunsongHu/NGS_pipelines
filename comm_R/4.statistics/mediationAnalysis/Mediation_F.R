######################################## 1. mediation analysis based on mediation package ########################################

#treat_key="Mother_OverweightObesity_VS_Normal"
#mediation="X52984"
#outcome="BirthBMI_age_z"

#rm(treat_key);rm(mediation);rm(outcome)

cal_mediation_one=function(indata,treat_key,mediation,mediation_type,outcome,outcome_type,cov_in){
  cov_all=c(treat_key,mediation,outcome,cov_in)
  indata=indata[apply(apply(indata[cov_all],2,function(x){ifelse(is.na(x),1,0)}),1,sum)<1,]
  n=nrow(indata)
  #for continuous outcome and continuous mediation
  if(outcome_type=="con" & mediation_type=="con"){
    indata1=data.frame(
      obs=1:nrow(indata),
      m=unlist(indata[mediation]),
      y=as.numeric(unlist(indata[outcome])),
      x=unlist(indata[treat_key]),
      stringsAsFactors = F
    )
    
    indata$obs=1:nrow(indata)
    data_cov=indata[c("obs",cov_in)]
    indata1=indata1 %>% left_join(data_cov)
    
    formula_m=paste0("m ~ x+",  paste0(cov_in,collapse = "+"))
    formula_y=paste0("y ~ x+m+",paste0(cov_in,collapse = "+"))
    
    model.m <- glm(formula_m, data = indata1)
    model.y <- glm(formula_y, data = indata1)
    
    set.seed(123456)
    fit_med_out=mediation::mediate(model.m, model.y, sims = sim_time, treat = "x",mediator = "m")
  }
  
  #for categorical outcome and continuous mediation
  if(outcome_type=="cat" & mediation_type=="con"){
    indata1=data.frame(
      obs=1:nrow(indata),
      m=unlist(indata[mediation]),
      y=as.numeric(factor(unlist(indata[outcome])))-1,
      x=unlist(indata[treat_key]),
      stringsAsFactors = F
    )
    
    indata$obs=1:nrow(indata)
    data_cov=indata[c("obs",cov_in)]
    indata1=indata1 %>% left_join(data_cov)
    
    formula_m=paste0("m ~ x+",  paste0(cov_in,collapse = "+"))
    formula_y=paste0("y ~ x+m+",paste0(cov_in,collapse = "+"))
    
    model.m <- glm(formula_m, data = indata1)
    model.y <- glm(formula_y, family = 'binomial',data = indata1)
    
    set.seed(123456)
    fit_med_out=mediation::mediate(model.m, model.y, sims = sim_time, treat = "x",mediator = "m")
  }
  
  #for categorical outcome and categorical mediation
  if(outcome_type=="cat" & mediation_type=="cat"){
    indata1=data.frame(
      obs=1:nrow(indata),
      m=as.numeric(factor(unlist(indata[mediation])))-1,
      y=as.numeric(factor(unlist(indata[outcome])))-1,
      x=unlist(indata[treat_key]),
      stringsAsFactors = F
    )
    
    indata$obs=1:nrow(indata)
    data_cov=indata[c("obs",cov_in)]
    indata1=indata1 %>% left_join(data_cov)
    
    
    formula_m=paste0("m ~ x+",  paste0(cov_in,collapse = "+"))
    formula_y=paste0("y ~ x+m+",paste0(cov_in,collapse = "+"))
    
    model.m <- glm(formula_m, family = 'binomial',data = indata1)
    model.y <- glm(formula_y, family = 'binomial',data = indata1)
    
    set.seed(123456)
    fit_med_out=mediation::mediate(model.m, model.y, sims = sim_time, treat = "x",mediator = "m")
  }
  
  #summarise output
  #suummary(fit_mediation)
  out=data.frame(
    outcome=outcome,
    Exposure=treat_key,
    Mediation=mediation,
    N_used=n,
    
    effect_total=paste0(sprintf("%.4f",round(fit_med_out$tau.coef,4)), " (",paste0(sprintf("%.4f",round(fit_med_out$tau.ci,4)),collapse = " - "),")"),
    p_total=fit_med_out$tau.p,
    
    effect_direct=paste0(sprintf("%.4f",round(fit_med_out$z.avg,4)), " (",paste0(sprintf("%.4f",round(fit_med_out$z.avg.ci,4)),collapse = " - "),")"),
    p_direct=fit_med_out$z.avg.p,
    
    effect_indirect=paste0(sprintf("%.4f",round(fit_med_out$d.avg,4)), " (",paste0(sprintf("%.4f",round(fit_med_out$d.avg.ci,4)),collapse = " - "),")"),
    p_indirect=fit_med_out$d.avg.p,
    
    mediation_proportion=paste0(
      sprintf("%.1f",round(fit_med_out$n.avg,4)*100),"%",
      " (",
      paste0(
        paste0(sprintf("%.1f",round(fit_med_out$n.avg.ci,4)*100),"%"),
        collapse = " - "),")")
    ,
    p_proportion=fit_med_out$n.avg.p,
    
    Cov=paste0(cov_in,collapse = "+"),
    stringsAsFactors = F
  )
  list(out=out,fit_med_out=fit_med_out)
  
}












































