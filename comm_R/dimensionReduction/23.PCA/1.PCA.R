

#PCA using data frame ----
run_pca_on_df=function(df_in,var_x_list){
  df_in_for_analysis=(df_in %>% mutate(obs=1:n()))[c("obs",var_x_list)]
  
  df_in_for_analysis$Missing=
    apply(df_in_for_analysis[2:ncol(df_in_for_analysis)],1,
          function(x){ifelse(any(is.na(x)),"HasMissing","NoMissing")})
  
  df_in_for_analysis1=df_in_for_analysis %>% filter(Missing=="NoMissing")
  
  fit_pca=prcomp(df_in_for_analysis1[var_x_list], center = TRUE,scale. = TRUE)
  
  print(summary(fit_pca))
  
  df_pc=as.data.frame(predict(fit_pca)) %>% mutate(obs=df_in_for_analysis1$obs)
  
  df_out=df_in_for_analysis %>% left_join(df_pc)
  df_out[c("obs","Missing",var_x_list)]=NULL
  
  return(list(fit_pca=fit_pca,df_pc=df_out))
}
