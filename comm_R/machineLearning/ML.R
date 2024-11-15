# run boruta ----------------------------------------
library(Boruta)

run_boruta=function(matrix_in,diag_in,max_run){
  data_with_diagnosis=as.data.frame(t(matrix_in)) %>%
    mutate(COH_sample=colnames(matrix_in)) %>% 
    left_join(diag_in[c("COH_sample","diag")]) %>%
    mutate(COH_sample=NULL,diag=as.factor(diag)) %>%  drop_na()
  
  print(dim(data_with_diagnosis))
  
  set.seed(10)
  boruta_out = Boruta(diag~., data = data_with_diagnosis, doTrace = 2,maxRuns=max_run)
  boruta_out_fixed = TentativeRoughFix(boruta_out)
  
  importance_boruta=attStats(boruta_out)
  importance_boruta_fixed=attStats(boruta_out_fixed)
  
  data_boruta_confirmed=matrix_in[row.names(importance_boruta)[importance_boruta$decision=="Confirmed"],]
  data_boruta_confirmed_fixed=matrix_in[row.names(importance_boruta_fixed)[importance_boruta_fixed$decision=="Confirmed"],]
  
  return(list(
    fit_boruta=boruta_out,
    importance_boruta=importance_boruta,importance_boruta_fixed=importance_boruta_fixed,
    data_boruta_confirmed=data_boruta_confirmed,data_boruta_confirmed_fixed=data_boruta_confirmed_fixed))
}

# KNN prediction -------------------------------------


KNN_pred_df=function(indata,outcome,varX_list,id_var,KNN_K){
  formula_in=formula(paste0(outcome,"~",paste0(varX_list,collapse = "+")))
  
  
  knnFit = caret::train(formula_in, data = indata, method = "knn",
                        preProcess = c("center", "scale"),
                        tuneGrid = expand.grid(k = c(KNN_K)))
  
  knn_pred = predict(knnFit) %>% as.vector()
  outdata=indata[c(id_var,outcome)] %>% mutate(diag_pred_knn=knn_pred)
  outdata
}

knn_pred=function(indata){
  KNN_K=3
  knnFit = train(diag ~X+Y, data = indata, method = "knn",
                 preProcess = c("center", "scale"),
                 tuneGrid = expand.grid(k = c(KNN_K)))
  
  knn_pred = predict(knnFit) %>% as.vector()
  outdata=indata %>% select(COH_sample,diag) %>% mutate(diag_pred_knn=knn_pred)
  outdata
}

knn_pred_file=function(file,diag_in){
  df_tsne=read_tsv(file) %>% mutate(group=paste0(gene_n,"_",perplexity)) %>% left_join(diag_in)
  
  indata=df_tsne
  groups=unique(indata$group)
  
  for(group in groups){
    g_index=match(group,groups)
    print(g_index)
    indata_one=indata[indata$group==group,]
    predict_one=knn_pred(indata_one)
    names(predict_one)[3]=paste0("pred_",group)
    if(g_index==1){predict_out=predict_one}
    if(g_index>1){predict_out=predict_out %>% left_join(predict_one)}
  }
  predict_out$knn_pred=apply(predict_out[3:ncol(predict_out)],1,function(x){names(sort(table(x),decreasing=TRUE)[1])})
  predict_out
}

knn_pred_forDf=function(indata,diag_in){
  indata=indata %>% mutate(group=paste0(gene_n,"_",perplexity)) %>% left_join(diag_in)
  groups=unique(indata$group)
  KNN_K=3
  for(group in groups){
    g_index=match(group,groups)
    print(g_index)
    indata_one=indata[indata$group==group,]
    predict_one=knn_pred(indata_one)
    names(predict_one)[3]=paste0("pred_",group)
    if(g_index==1){predict_out=predict_one}
    if(g_index>1){predict_out=predict_out %>% left_join(predict_one)}
  }
  predict_out
}


cal_classifcation_rate=function(indata,var_y,var_predicted){
  data_freq_classcification=as.data.frame(table(indata[var_y]==unlist(indata[var_predicted]),unlist(indata[var_y])))
  data_freq_classcification1=tidyr::spread(data_freq_classcification,Var1,Freq) 
  names(data_freq_classcification1)[1:3]=c(var_y,"False","True")
  data_freq_classcification2=data_freq_classcification1%>%
    mutate(
      N_subtype=False +True,
      Right=paste0(True," (",sprintf("%.2f",100*True/N_subtype),"%)"),
      Wrong=paste0(False," (",sprintf("%.2f",100*False/N_subtype),"%)"),
      True=NULL,False=NULL
    )
  data_freq_classcification2
}


















