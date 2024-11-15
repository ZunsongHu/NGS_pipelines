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

########################################  gee, one x with 2 outcomes at the same time ########################################


#names(data_analysis)[1:20]


#indata=data_analysis
#var_id="studyid"
#var_outcome1="Mother_OverweightObesity_VS_Normal"
#var_outcome2="Child_cv4_OverweightObesity_VS_Noraml"
#var_x_one="X38768"
#adj=c("mother_age","mother_insurance","mother_edu","mother_tobacco","mother_alcohol","parity")
#digit_n=2  

#var_x_list=meta_name[1:10]

#x=cal_gee_OneVarX_TwoVarOutcome(data_analysis,"studyid","Mother_OverweightObesity_VS_Normal","Child_cv4_OverweightObesity_VS_Noraml",meta_name[1:10],adj,2)
  

cal_gee_OneVarX_TwoVarOutcome=function(indata,var_id,var_outcome1,var_outcome2,var_x_list,adj,digit_n){
  bind_rows(lapply(var_x_list, function(var_x_one){
    if(is.na(adj)[1]|adj[1]==""){indata=data_analysis[c(var_id,var_outcome1,var_outcome2,var_x_one)]}
    if(!(is.na(adj)[1]|adj[1]=="")){indata=data_analysis[c(var_id,var_outcome1,var_outcome2,var_x_one,adj)]}
    
    indata=indata %>% na.omit()
    indata$id=unlist(indata[var_id])
    
    if(!(is.na(adj)[1]|adj[1]=="")){
      fit_formula_half=formula(paste0(var_x_one," ~ ",paste0(c(adj),collapse = "+")))
      fit_formula_full=formula(paste0(var_x_one," ~ ",paste0(c(var_outcome1,var_outcome2,adj),collapse = "+")))
    }
    
    fit_gee1=geepack::geeglm(fit_formula_half, data = indata,id = id,corstr = "exchangeable")
    fit_gee2=geepack::geeglm(fit_formula_full , data = indata,id = id,corstr = "exchangeable")
    fit_anova=anova(fit_gee1,fit_gee2)
    
    style_number=paste0("%.",digit_n,"f")
    
    out_one=data.frame(
      Outcome1=var_outcome1,
      Outcome2=var_outcome2,
      variable=var_x_one,
      Adj=paste0(c(adj),collapse = "+"),
      
      N_used_model1=length(fit_gee1$weights),
      Formula_model1=paste0(as.character(fit_gee1$formula)[2:3],collapse = " ~ "),
      Beta_model1=paste0(
        paste0(names(fit_gee1$coefficients),": ",
               sprintf(style_number,fit_gee1$coefficients)
        ),collapse = ";"),
      
      N_used_model2=length(fit_gee2$weights),
      Formula_model2=paste0(as.character(fit_gee2$formula)[2:3],collapse = " ~ "),
      Beta_model2=paste0(
        paste0(names(fit_gee2$coefficients),": ",
               sprintf(style_number,fit_gee1$coefficients)
        ),collapse = ";"),
      
      
      anova_chisq=fit_anova$X2,
      anova_df1_p_raw=fit_anova$`P(>|Chi|)`,
      anova_df1_p=convert_p(fit_anova$`P(>|Chi|)`),
      anova_df2_p_raw=pchisq(fit_anova$X2, df=2, lower.tail=FALSE),
      anova_df2_p=convert_p(pchisq(fit_anova$X2, df=2, lower.tail=FALSE)),
      stringsAsFactors = F
    )
    out_one
    
  }))
}







