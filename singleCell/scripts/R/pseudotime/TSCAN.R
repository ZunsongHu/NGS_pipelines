#Pseudotime analysis  ------------------------------------------------------------------------------------------------------------------------------------------------------

get_pseudotime=function(obj_in,var_cluster,pseduotime_start,variable_n=2000){
  library(scater)
  library(TSCAN)
  
  #data prepare
  df_info=obj_in@meta.data
  
  df_info$group_mst=unlist(df_info[var_cluster])
  
  df_info=df_info %>% arrange(group_mst)
  
  
  matrix_cont=GetAssayData(obj_in,assay = "RNA",slot="counts")
  matrix_cont=matrix_cont[,row.names(df_info)]
  
  #rerun normalization to pca
  message("Rerun normalization to pca...\n")
  obj_1=CreateSeuratObject(matrix_cont,meta.data = df_info)
  obj_1=Seurat::NormalizeData(obj_1,verbose =F)
  obj_1=Seurat::FindVariableFeatures(obj_1,nfeatures=variable_n,verbose =F)
  obj_1=Seurat::ScaleData(obj_1,verbose =F)
  obj_1=RunPCA(obj_1,verbose =F)
  # obj_1=RunTSNE(obj_1,dims = 1:50,verbose =F)
  Idents(obj_1)=var_cluster
  levels(obj_1)=sort(unlist(unique(obj_1[[var_cluster]])))
  
  #get as.SingleCellExperiment obj
  obj_in1=as.SingleCellExperiment(obj_1)
  
  #get the cluster centroids by averaging the coordinates of its member cells
  centroids=reducedDim(aggregateAcrossCells(obj_in1,ids=obj_in1[[var_cluster]]),"PCA")
  
  #form the minimum spanning tree across those centroids
  mst=TSCAN::createClusterMST(centroids,cluster=NULL)
  #obtain a pseudotime ordering by projecting the cells onto the MST. That is move each cell onto the most closet edge of the MST
  map.tscan=mapCellsToEdges(obj_in1,mst=mst,clusters=obj_in1[[var_cluster]],use.dimred="PCA")
  #the pseduotime is then calculated as the distance along the MST to this new position from a root node with orderCells(), Need to pick a root node.
  tscan.pseudo=orderCells(map.tscan,mst,start = pseduotime_start)
  #the common pseudo
  common.pseudo=averagePseudotime(tscan.pseudo)
  
  #DE genes along the trajectory
  DE_pseudotime=as.data.frame(testPseudotime(obj_in1,pseudotime=common.pseudo))
  
  DE_pseudotime$feature=row.names(DE_pseudotime) 
  DE_pseudotime=DE_pseudotime%>% arrange(desc(abs(logFC)))
  
  list(tscan.pseudo=tscan.pseudo,common.pseudo=common.pseudo,DE_pseudotime=DE_pseudotime)
}

