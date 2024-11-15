######################################## Extract results from meta objects ########################################

meta_extract=function(meta_object,digit_n){
  one=data.frame(study_number=meta_object$k,
                 tau=round(meta_object$tau,4),
                 I2=paste0(round(meta_object$I2,4)*100,"%"),
                 Q=round(meta_object$Q,4),Q_test_p=signif(2*(1-pchisq(meta_object$Q,meta_object$df.Q))),
                 
                 effect_fixed=round(meta_object$TE.fixed,digit_n),
                 effect_se_fixed=round(meta_object$seTE.fixed,digit_n),
                 effect_ll_fixed=round(meta_object$lower.fixed,digit_n),
                 effect_ul_fixed=round(meta_object$upper.fixed,digit_n),
                 fixed_z=round(meta_object$zval.fixed,4),
                 fixed_p=signif(meta_object$pval.fixed,digits = 3),
                 
                 effect_random=round(meta_object$TE.random,digit_n),
                 effect_se_random=round(meta_object$seTE.random,digit_n),
                 effect_ll_random=round(meta_object$lower.random,digit_n),
                 effect_ul_random=round(meta_object$upper.random,digit_n),
                 random_z=round(meta_object$zval.random,4),
                 random_p=signif(meta_object$pval.random,digits = 3),
                 
                 egger_test_p=NA,begg_test_p=NA,
                 stringsAsFactors = F
  )
  
  if (meta_object$k>=4){
    egger_test_p=signif(meta::metabias(meta_object,method.bias = c("linreg"),k.min = 4)$p.value,4)
    begg_test_p=signif(meta::metabias(meta_object,method.bias = c("rank"),k.min = 4)$p.value,4)
    one[(length(one)-1)]=egger_test_p
    one[length(one)]=begg_test_p
  }
  one
}




######################################## meta analysis beta se ########################################


cal_meta_betaSe_one=function(beta_values,se_values){
  meta_fit = meta::metagen(beta_values,se_values,comb.fixed = TRUE,comb.random = T)
  meta_extract(meta_fit,6)
}


######################################## meta analysis or ########################################


#indata=data_biomarker
#var_beta="lnrr"
#var_se="selnrr"
#var_label="author_year"
#var_subgroup="index"

#dir_forest_plot="2.data_temp/data_for_stata1/forest_plot/"

meta_beta_se=function(indata,var_beta,var_se,var_label,var_subgroup,var_type,filelable){
  value_beta=unlist(indata[var_beta])
  value_se=unlist(indata[var_se])
  value_label=unlist(indata[var_label])
  value_subgroup=unlist(indata[var_subgroup])
  
  meta_object=meta::metagen(value_beta,value_se,comb.fixed = TRUE,comb.random = T,sm = var_type,
                            studlab=value_label,
                            byvar = value_subgroup)
  
  pdf(paste0(dir_forest_plot,filelable,"_",var_subgroup,".pdf"),width=8,height =16)
  
  meta::forest(meta_object, layout="JAMA",
               overall=F,overall.hetstat=F,print.tau2=F,print.pval.Q=F,
               col.study="black",col.square="gray",
               col.diamond.fixed="blue",col.diamond.random="orange",
               col.by="black")
  dev.off()
  
}





















