# functions ------------------------------------------------------------------------------------------- ----

get_celltype_df=function(inObj,df_ref){
  for_singleR=GetAssayData(inObj)
  
  pred_celltype <- SingleR::SingleR(test = for_singleR, ref = df_ref, assay.type.test=1,labels = df_ref$cellType)
  df_label_pred=data.frame(
    barcode=row.names(inObj@meta.data),
    pred_celltype=pred_celltype$labels,
    stringsAsFactors = F)
  df_label_pred
}

#SingleR Annotation -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# inObj=obj_QC
# ref_singleR=ref_BALL
# var_group="diag"

create_singleR_ref=function(exp,pheno,var_group){
  pheno1=pheno[var_group]
  
  if(!all(colnames(exp)==row.names(pheno1))){stop("Exp colname not equal to row names of Pheno")}
  
  objOut = SummarizedExperiment(list(logcounts=exp), colData = pheno1)
  objOut
}



# inObj=obj_QC
# ref_singleR=ref_BALL
# var_group="diag"

run_singleR=function(inObj,ref_singleR,var_group,
                     transform="SCT",method="classic"){
  #using top variable features
  # top_variable_features=head(VariableFeatures(inObj), top_variable_features_n)
  # for_singleR=GetAssayData(inObj)
  # for_singleR=for_singleR[top_variable_features,]
  
  if(transform=="SCT"){
    #using all the features
    for_singleR=GetAssayData(inObj[["SCT"]], slot = "data")
  }
  
  if(transform=="ordinary"){
    #using all the features
    for_singleR=GetAssayData(inObj)
  }
  
  print(dim(for_singleR))
  
  value_group=unlist(ref_singleR@colData[var_group])
  
  #cell level annotation
  print("cell level annotation")
  out_singleR = SingleR::SingleR(test = for_singleR, ref = ref_singleR,labels = value_group)
  
  #cluster level annotation
  print("cluster level annotation")
  out_singleR_cluster = SingleR::SingleR(test = for_singleR, ref = ref_singleR,labels = value_group,clusters=inObj$seurat_clusters)
  
  #get output 
  df_singleR_out=data.frame(
    barcode=row.names(inObj@meta.data),
    seurat_clusters=inObj$seurat_clusters,
    labels=out_singleR$labels,
    pruned.labels=out_singleR$pruned.labels,
    stringsAsFactors = F) %>%
    left_join(data.frame(
      seurat_clusters=sort(unique(inObj$seurat_clusters)),
      labels_cluster=out_singleR_cluster$labels,
      pruned.labels_cluster=out_singleR_cluster$pruned.labels,
      stringsAsFactors = F
    )
    )
  
  return(list(out_singleR=out_singleR,out_singleR_cluster=out_singleR_cluster,df_singleR_out=df_singleR_out))
}

#' run_singleR_clusterOnly
#'
#' @param inObj 
#' @param ref_singleR 
#' @param var_group 
#' @param transform 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
run_singleR_clusterOnly=function(inObj,ref_singleR,var_group,transform="ordinary",var_cluster="seurat_clusters"){
  #using top variable features
  # top_variable_features=head(VariableFeatures(inObj), top_variable_features_n)
  # for_singleR=GetAssayData(inObj)
  # for_singleR=for_singleR[top_variable_features,]
  
  if(transform=="SCT"){
    #using all the features
    for_singleR=GetAssayData(inObj[["SCT"]], slot = "data")
  }
  
  if(transform=="ordinary"){
    #using all the features
    for_singleR=GetAssayData(inObj)
  }
  
  print(dim(for_singleR))
  
  value_group=unlist(ref_singleR@colData[var_group])
  
  #cluster level annotation
  print("cluster level annotation")
  value_cluster=unlist(inObj@meta.data[var_cluster])
  
  out_singleR_cluster = SingleR::SingleR(test = for_singleR, ref = ref_singleR,labels = value_group,clusters=value_cluster)
  
  #get output 
  df_singleR_out=data.frame(
    barcode=row.names(inObj@meta.data),
    var_cluster=value_cluster,
    stringsAsFactors = F) %>%
    left_join(data.frame(
      var_cluster=sort(unique(value_cluster)),
      labels_cluster=out_singleR_cluster$labels,
      pruned.labels_cluster=out_singleR_cluster$pruned.labels,
      stringsAsFactors = F
    )
    )
  
  return(list(out_singleR_cluster=out_singleR_cluster,df_singleR_out=df_singleR_out))
}




# inObj=obj_QC
# ref_singleR1=ref_HPCA
# ref_singleR2=ref_BlineageRef_30types
# var_group1="label.final"
# var_group2="celltype"


run_singleR_doubleRef=function(inObj,ref_singleR1,ref_singleR2,var_group1,var_group2){
  #using all the features
  # for_singleR=GetAssayData(inObj[["SCT"]], slot = "data")
  for_singleR=GetAssayData(inObj)
  
  print(dim(for_singleR))
  
  value_group1=unlist(ref_singleR1@colData[var_group1])
  value_group2=unlist(ref_singleR2@colData[var_group2])
  
  #cell level annotation
  print("cell level annotation")
  out_singleR = SingleR::SingleR(test = for_singleR, 
                                 ref = list(ref1=ref_singleR1,ref2=ref_singleR2),
                                 labels = list(value_group1,value_group2))
  
  #cluster level annotation
  print("cluster level annotation")
  out_singleR_cluster = SingleR::SingleR(test = for_singleR, 
                                         ref = list(ref1=ref_singleR1,ref2=ref_singleR2),
                                         labels = list(value_group1,value_group2),
                                         clusters=inObj$seurat_clusters)
  
  #get output 
  df_singleR_out=data.frame(
    barcode=row.names(inObj@meta.data),
    seurat_clusters=inObj$seurat_clusters,
    labels=out_singleR$labels,
    pruned.labels=out_singleR$pruned.labels,
    stringsAsFactors = F) %>%
    left_join(data.frame(
      seurat_clusters=sort(unique(inObj$seurat_clusters)),
      labels_cluster=out_singleR_cluster$labels,
      pruned.labels_cluster=out_singleR_cluster$pruned.labels,
      stringsAsFactors = F
    )
    )
  return(list(out_singleR=out_singleR,out_singleR_cluster=out_singleR_cluster,df_singleR_out=df_singleR_out))
}


draw_annotation_plot=function(objIn,singleR_obj,label,cols_in,dir_out){
  
  objIn[[label]]=singleR_obj$pruned.labels
  
  DimPlot(objIn,group.by = label,label=T,repel = T,
          cols = cols_in[names(cols_in) %in% unlist(objIn@meta.data[label])])
  ggsave(paste0(dir_out,label,".png"),width=11,height = 10)
  ggsave(paste0(dir_out,label,".pdf"),width=11,height = 10)
  
  objIn
}

get_singleR_labels=function(in_singleR_obj){
  data.frame(
    barcode=in_singleR_obj@rownames,
    labels=in_singleR_obj$labels,
    pruned.labels=in_singleR_obj$pruned.labels,
    deltaFromMedian=SingleR::getDeltaFromMedian(in_singleR_obj),
    stringsAsFactors = F
  )
}

get_singleR_scoreMatrix=function(in_singleR_obj){
  matrix_score=bind_cols(
    data.frame(
      barcode=in_singleR_obj@rownames,
      labels=in_singleR_obj$labels,
      stringsAsFactors = F
    ),
    as.data.frame(in_singleR_obj$scores)
  )
  row.names(matrix_score)=matrix_score$barcode
  matrix_score
}

get_singleR_heatmapScore=function(in_singleR_obj){
  df_label=get_singleR_labels(in_singleR_obj) %>% 
    left_join(get_singleR_topCells(in_singleR_obj)) %>% 
    filter(labels==labelsFromMaxScore) %>% mutate(n=n()) %>% 
    group_by(labels) %>% mutate(n_g=n(),diag_per=round(n_g/n,4)) %>% 
    arrange(desc(n_g),desc(deltaFromMedian)) %>% ungroup() %>% mutate(obs=1:n())
  
  df_label_unique=df_label %>% dplyr::select(labels,n_g,diag_per) %>% arrange(desc(n_g)) %>% distinct() %>% 
    mutate(label_=paste0(letters[1:n()],".",labels),
           per=paste0(sprintf("%.1f",diag_per*100),"%"),
           label_per=paste0(labels," (",per,")"))
  
  names(df_label_unique)[1]="variable"
  
  matrix_score_=df_label %>% 
    left_join(get_singleR_scoreMatrix(in_singleR_obj) %>% mutate(labels=NULL)) %>% 
    arrange(obs)
  
  scale_0to1=function(x){
    x=as.numeric(x)
    (x-min(x))/(max(x)-min(x))}
  
  df_score_=matrix_score_[c("barcode","obs","labels",df_label_unique$variable)] %>% 
    reshape2::melt(id.vars=c("barcode","obs","labels")) %>% left_join(df_label_unique) %>% 
    group_by(barcode) %>% 
    mutate(value1=scale_0to1(value)) %>% ungroup()
  
  df_score_
}

get_singleR_topCells=function(in_singleR_obj,topPercentileOfDeltaToKeep=0.9,remove_valueLabelNonConsistent=F){
  #get label using max score
  df_score=as.data.frame(in_singleR_obj$scores)
  df_score$barcode=in_singleR_obj@rownames
  df_score1=df_score %>% reshape2::melt(id.vars=c("barcode")) %>% 
    group_by(barcode) %>% filter(value==max(value))
  names(df_score1)[2:3]=c("labelsFromMaxScore","maxScore")
  
  #get top cells using delta
  cutoff_delta_precentile=1-topPercentileOfDeltaToKeep
  df_delta=get_singleR_labels(in_singleR_obj) %>% 
    group_by(labels) %>% 
    mutate(
      delta_cutoff=quantile(deltaFromMedian,cutoff_delta_precentile),
      label_top=ifelse(deltaFromMedian>delta_cutoff,as.character(labels),NA)
    ) %>% ungroup()
  
  #get output 
  df_out=df_score1 %>% left_join(df_delta %>% dplyr::select(barcode,label_top))
  
  if(remove_valueLabelNonConsistent){
    df_out=df_out %>% 
      mutate(label_top=ifelse(label_top==labelsFromMaxScore,label_top,NA))
  }
  
  df_out
}


draw_barplot=function(x,file_out_prefix,width=5,height=10){
  x[is.na(x)]="NA";  n=length(x)
  
  df_freq=data.frame(table(x)) %>% arrange(desc(Freq)) %>%
    mutate(Group=as.character(x),
           ratio=Freq/n,
           Percentage=paste0(Group,": ",Freq," (",sprintf("%.2f",ratio*100),"% )"))
  
  ggplot(df_freq,aes(x=ratio,y=Percentage)) + geom_bar(stat="identity") +
    theme_classic()
  ggsave(paste0(file_out_prefix,".png"),width=width,height = height)
  ggsave(paste0(file_out_prefix,".pdf"),width=width,height = height)
}

























