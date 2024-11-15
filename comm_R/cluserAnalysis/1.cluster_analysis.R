
#################### cal_best_cluster_n ######################



cal_best_cluster_n=function(indata,var_in){
  indata1=na.omit(indata[var_in])
  
  print("start hierarchical")
  set.seed(123)
  nbclust_fit_hierarchical=NbClust::NbClust(indata1, distance = "euclidean",
                                            min.nc = 2, max.nc = 10, 
                                            method = "complete", index ="all")
  
  nbclust_fit_hierarchical_best_n_data=nbclust_fit_hierarchical$Best.nc
  nbclust_fit_hierarchical_best_n=as.numeric(as.character(data.frame(table(nbclust_fit_hierarchical_best_n_data[1,]))[,1]
                                                          [which.max(data.frame(table(nbclust_fit_hierarchical_best_n_data[1,]))[,2])]))
  
  print("start kmeans")
  set.seed(123)
  nbclust_fit_kmeans=NbClust::NbClust(indata1, distance = "euclidean",
                                      min.nc = 2, max.nc = 10, 
                                      method = "kmeans", index ="all")
  
  nbclust_fit_kmeans_best_n_data=nbclust_fit_kmeans$Best.nc
  nbclust_fit_kmeans_best_n=as.numeric(as.character(data.frame(table(nbclust_fit_kmeans_best_n_data[1,]))[,1]
                                                    [which.max(data.frame(table(nbclust_fit_kmeans_best_n_data[1,]))[,2])]))
  
  
  #factoextra::fviz_nbclust(nbclust_fit_hierarchical) 
  #factoextra::fviz_nbclust(nbclust_fit_kmeans) 
  
  data_n=data.frame(
    var_in=paste0(var_in,collapse = "+"),
    hierarchical_best_n=nbclust_fit_hierarchical_best_n,
    kmeans_best_n=nbclust_fit_kmeans_best_n,
    stringsAsFactors = F
  )
  
  data_n
}



#################### k-means ######################


cal_cluster_kmeans=function(indata,var_in,group,var_id,var_cluster){
  indata1=na.omit(indata[c(var_id,var_in)])
  set.seed(123)
  fit_kemans=kmeans(indata1[2:ncol(indata1)], centers = group, nstart = 25)
  
  pdf(paste0(plot_dir,var_cluster,".pdf"))
  print(factoextra::fviz_cluster(fit_kemans,data=indata1[2:ncol(indata1)]))
  dev.off()
  
  data_out=data.frame(
    id=unlist(indata1[,1]),
    cluster=fit_kemans$cluster,
    stringsAsFactors = F
  )
  names(data_out)=c(var_id,var_cluster)
  data_out
}

#################### hierarchical ######################



cal_cluster_hierarchical=function(indata,var_in,group,var_id,var_cluster){
  
  indata1=na.omit(indata[c(var_id,var_in)])
  
  d= dist(indata1[2:ncol(indata1)], method = "euclidean")
  set.seed(123)
  hc1 = hclust(d, method = "complete" )
  
  pdf(paste0(plot_dir,var_cluster,"_tree_.pdf"))
  plot(hc1,cex=0.6)
  rect.hclust(hc1, k = group, border = 2:(group+1))
  dev.off()
  
  sub_grp = cutree(hc1, k = group)
  
  pdf(paste0(plot_dir,var_cluster,"_.pdf"))
  print(factoextra::fviz_cluster(list(data = indata1[2:ncol(indata1)], cluster = sub_grp)))
  dev.off()
  
  data_out=data.frame(
    id=unlist(indata1[,1]),
    cluster=sub_grp,
    stringsAsFactors = F
  )
  names(data_out)=c(var_id,var_cluster)
  data_out
}







########### get cluster groups elbew method
#set.seed(123)
#wss_fit=fviz_nbclust(indata1, FUN = hcut, method = "wss")
#wss_fit_data=wss_fit$data
#wss_fit_n=as.numeric(wss_fit_data$clusters[which.max(wss_fit_data$y)])

#set.seed(123)
#silhouette_fit=fviz_nbclust(indata1, FUN = hcut, method = "silhouette")
#silhouette_fit_data=silhouette_fit$data
#silhouette_fit_n=as.numeric(silhouette_fit_data$clusters[which.max(silhouette_fit_data$y)])

#set.seed(123)
#gap_fit=fviz_nbclust(indata1, hcut, nstart = 25,  method = "gap_stat", nboot = 50)+
#  labs(subtitle = "Gap statistic method")
#gap_fit_data=gap_fit$data
#gap_fit_n=as.numeric(gap_fit_data$clusters[which.max(gap_fit_data$gap)])


#fit_clValid= clValid::clValid(as.matrix(indata1), nClust = 2:10, 
#                              clMethods = c("hierarchical","kmeans","pam"), validation = "internal")


#fit_clValid_summary<-clValid::summary(fit_clValid)
#clValid::plot(fit_clValid)