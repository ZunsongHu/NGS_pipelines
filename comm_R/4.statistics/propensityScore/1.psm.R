######################################## calculate psm ##################################

cal_psm=function(indata,treatment,var_key){
  indata$outcome=factor(unlist(indata[treatment]))
  indata$outcome=relevel(indata$outcome,ref =levels(indata$outcome)[1])
  formula_fit=formula(paste0('outcome',"~",paste0(var_key,collapse = "+")))
  fit_logit=glm(formula_fit,family = "binomial",data=indata)
  predict(fit_logit,type="response")
}

######################################## generate psm matched data ##################################

generate_psm_matched_data=function(indata,treatment,var_key,ratio){
  indata$outcome=as.numeric(factor(unlist(indata[treatment])))-1
  if(table(indata$outcome)[1]< table(indata$outcome)[2]){indata$outcome=ifelse(indata$outcome==0,1,0)}
  
  formula_fit=formula(paste0('outcome',"~",paste0(var_key,collapse = "+")))
  mod_match = MatchIt::matchit(formula_fit,method = "nearest", data = indata,ratio=ratio)
  
  data_matched=MatchIt::match.data(mod_match)
  data_matched$outcome=NULL
  data_matched
}
