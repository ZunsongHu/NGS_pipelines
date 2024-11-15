######################################## LCA one function ########################################
library(poLCA)

lca_one=function(indata,feature_in,group,plot_dir,plot_name,seed_one){
  x=complete.cases(indata[feature_in])
  indata=indata[x,]
  indata[feature_in]=indata[feature_in]+1
  
  f1=formula(paste0("cbind(",paste0(feature_in,collapse = ","),")~1"))
  
  
  #pdf(paste0(plot_dir,plot_name,".pdf"),width =20,height = 20)
  set.seed(3)
  par(mar=c(5.1, 4.1, 4.1, 10))
  fit_lca1=poLCA::poLCA(f1,indata,nclass=group,graphs = T,na.rm = T,maxiter = 5000,verbose=F)
  #dev.off()
  
  
  statistic=data.frame(
    label=paste0(plot_name),
    Group=group,
    n_used=fit_lca1$Nobs,
    iteration_n_at_converge=fit_lca1$numiter,
    Group_percentage=paste0(sprintf("%0.3f",round(fit_lca1$P,3)),collapse = "/"),
    bic=fit_lca1$bic,
    aic=fit_lca1$aic,
    stringsAsFactors = F)
  group=data.frame(fit_lca1$predclass,stringsAsFactors = F)
  names(group)=plot_name
  
  return(list(
    group=group,
    statistic=statistic
  ))
}

