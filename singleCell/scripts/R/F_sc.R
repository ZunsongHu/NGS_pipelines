options(warn=-1)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(SummarizedExperiment)
  library(stringr)
  library(scCustomize)
  
  # library(magrittr)
  # library(viridis)
  # library(qs)
})
options(warn=0)



#args <- commandArgs(trailingOnly = TRUE)
#print(args)

#file_tsv=args[1]

add_percentage=function(inObj,gene_list,label){
  gene_list_all=rownames(x = inObj)
  inObj[[label]] <- PercentageFeatureSet(inObj, features = gene_list_all[gene_list_all %in% gene_list])
  inObj
}

#FilterGenes -------------------------------------
FilterGenes <-
  function (object, min.value=1, min.cells = 0, genes = NULL) {
    genes.use <- rownames(object)
    if (!is.null(genes)) {
      genes.use <- intersect(genes.use, genes)
      object@data <- GetAssayData(object)[genes.use, ]
    } else if (min.cells > 0) {
      num.cells <- Matrix::rowSums(GetAssayData(object) > min.value)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      object = object[genes.use, ]
    }
    object <- LogSeuratCommand(object = object)
    return(object)
  }




#create_seurat_object------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# min.cells, Include features detected in at least this many cells. Will subset the counts matrix as well. 
# min.features, Include cells where at least this many features are detected.
create_seurat_object=function(data_dir,
                              project_name="Seurat_obj",
                              min.cells=0,
                              min.features=0,
                              gene.column=2){
  print("reading data")
  data = Read10X(data.dir = data_dir,gene.column=gene.column)
  print("create object")
  seurat_object=CreateSeuratObject(counts = data, project = project_name, 
                                   min.cells = min.cells, min.features = min.features)
  return(seurat_object)
}


#Basic -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
get_barcode_df=function(object_in){
  data.frame(barcode=row.names(object_in@meta.data),stringsAsFactors = F)
}

get_embedding_df=function(inobj,reduction_method){
  df_embed = as.data.frame(Embeddings(inobj[[reduction_method]]))
  df_embed$barcode=row.names(df_embed)
  df_embed
}

get_MatrixCount_from_obj=function(inObj,outfile=NULL){
  matrix_count=as.matrix(GetAssayData(inObj[['RNA']],slot = "counts"))
  
  matrix_count_t=as.data.frame(as.matrix(t(matrix_count)))
  
  df_rownames=data.frame(cells=row.names(matrix_count_t),stringsAsFactors = F)
  
  matrix_out=bind_cols(df_rownames,matrix_count_t)
  
  if(!is.null(outfile)){
    write.csv(matrix_out,outfile,quote = F,row.names = F)
  }
  matrix_out
}

obj_subset_random=function(obj,n=2000){
  sampled_cells =sample(Cells(obj),n)
  obj_ = subset(obj_, cells = sampled_cells)
  obj_
}


#run sc analysis one sample -----------------
run_sc_analysis_toumap=function(inObj,label,MTgene,ribosomalGene,housekeepinngGene,dims){
  #add qc measurement
  inObj=add_percentage(inObj,MTgene,"percent.mt")
  inObj=add_percentage(inObj,ribosomalGene,"percent.ribosomal")
  inObj=add_percentage(inObj,housekeepinngGene,"percent.houseKeeping")
  
  #run pca
  inObj <- NormalizeData(inObj)
  inObj <- FindVariableFeatures(inObj, selection.method = "vst", nfeatures = 1000)
  inObj <- ScaleData(inObj)
  inObj <- RunPCA(inObj)
  inObj = CellCycleScoring(object = inObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
  
  #umap based on pca
  inObj <- RunUMAP(inObj, dims = dims)
  
  #run cluster
  inObj <- FindNeighbors(inObj, dims = dims)
  inObj <- FindClusters(inObj, resolution = 0.8)
  
  Idents(inObj)=label
  
  inObj
}


run_sc_analysis=function(inObj,label,MTgene,ribosomalGene,housekeepinngGene,dims){
  #add qc measurement
  inObj=add_percentage(inObj,MTgene,"percent.mt")
  inObj=add_percentage(inObj,ribosomalGene,"percent.ribosomal")
  inObj=add_percentage(inObj,housekeepinngGene,"percent.houseKeeping")
  
  #run pca
  inObj <- NormalizeData(inObj)
  inObj <- FindVariableFeatures(inObj, selection.method = "vst", nfeatures = 1000)
  inObj <- ScaleData(inObj)
  inObj <- RunPCA(inObj)
  inObj = CellCycleScoring(object = inObj, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
  
  #umap based on pca
  inObj <- RunUMAP(inObj, dims = dims)
  
  #run cluster
  inObj <- FindNeighbors(inObj, dims = dims)
  inObj <- FindClusters(inObj, resolution = 0.8)
  
  #get singlet doublet --------
  annotations=inObj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)  
  nExp_poi <- round(0.075*nrow(inObj@meta.data))  
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  inObj <- doubletFinder_v3(inObj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj,  sct = FALSE)
  
  Singlet_Doublet=names(inObj@meta.data)[grepl("DF.classifications",names(inObj@meta.data))]
  Singlet_Doublet_score=names(inObj@meta.data)[grepl("pANN",names(inObj@meta.data))]
  
  inObj[["Singlet_Doublet"]]=unlist(inObj@meta.data[Singlet_Doublet])
  inObj[["Singlet_Doublet_score"]]=unlist(inObj@meta.data[Singlet_Doublet_score])
  
  #umap based on sctransform
  # inObj_st <- SCTransform(inObj)
  # inObj_st <- RunPCA(inObj_st)
  # inObj_st <- RunUMAP(inObj_st, dims = dims)
  
  Idents(inObj)=label
  
  inObj
}


#run tsne for one sample -----------------
run_tsne=function(in_object,dir_plot,method_ForVariableGene="dispersion",N_variable=500,N_pc=40){
  if(!dir.exists(dir_plot)){dir.create(dir_plot)}
  
  #Normalizing the data
  in_object <- NormalizeData(in_object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  #Identification of highly variable features (feature selection)
  in_object <- FindVariableFeatures(in_object, selection.method = method_ForVariableGene, nfeatures = N_variable)
  
  top_10 <- head(VariableFeatures(in_object), 10)
  
  plot1 <- VariableFeaturePlot(in_object)
  plot2 <- LabelPoints(plot = plot1, points = top_10, repel = TRUE)
  plot1 + plot2
  ggsave(file.path(dir_plot,"VariableFeaturePlot.png"),width = 20,height = 8)
  ggsave(file.path(dir_plot,"VariableFeaturePlot.pdf"),width = 20,height = 8)
  
  #Scaling the data
  all.genes <- head(VariableFeatures(in_object), N_variable)
  
  in_object <- ScaleData(in_object, features = all.genes)
  
  #Perform linear dimensional reduction PCA
  set.seed(10)
  in_object <- RunPCA(in_object, features = all.genes,npcs = 50)
  
  DimPlot(in_object, reduction = "pca",dims = c(1,2),pt.size=2)
  ggsave(file.path(dir_plot,"DimPlot_PCA.png"),width = 10,height = 8)
  ggsave(file.path(dir_plot,"DimPlot_PCA.pdf"),width = 10,height = 8)
  
  png(file.path(dir_plot,"DimHeatmap_PCA.png"),width = 4000,height = 3000,res = 300)
  par(oma=c(1,1,1,10))
  DimHeatmap(in_object, dims = 1:9, cells = 500, balanced = TRUE)
  dev.off()
  
  pdf(file.path(dir_plot,"DimHeatmap_PCA.pdf"),width = 12,height = 9)
  par(oma=c(1,1,1,10))
  DimHeatmap(in_object, dims = 1:9, cells = 500, balanced = TRUE)
  dev.off()
  
  ElbowPlot(in_object,ndims = 50)
  ggsave(file.path(dir_plot,"ElbowPlot_PCA.png"),width = 10,height = 8)
  ggsave(file.path(dir_plot,"ElbowPlot_PCA.pdf"),width = 10,height = 8)
  
  #Run non-linear dimensional reduction (UMAP/tSNE)
  set.seed(10)
  in_object <- RunUMAP(in_object, reduction = "pca", dims = 1:N_pc)
  set.seed(10)
  in_object <- RunTSNE(in_object, reduction = "pca", dims = 1:N_pc)
  
  #cell-type annotation
  for_singleR=GetAssayData(in_object)
  
  pred_celltype <- SingleR(test = for_singleR, ref = hpca.se, assay.type.test=1,
                           labels = hpca.se$label.fine)
  
  df_label_pred=data.frame(pred_celltype=pred_celltype$labels,stringsAsFactors = F)
  df_label_pred = df_label_pred %>% mutate(
    pred_cell1=word(pred_celltype,1,sep=":"),
    pred_celltype_new=pred_celltype
  )
  #pred_celltype_new=ifelse(pred_cell1=="B_cell",pred_celltype,"Other")
  
  in_object[["SingleR.labels"]] <- df_label_pred$pred_celltype_new
  
  Idents(in_object)=df_label_pred$pred_celltype_new
  
  in_object1 <- subset(in_object, subset = SingleR.labels !="Other")
  
  DimPlot(in_object, reduction = "umap",pt.size = 1)
  ggsave(file.path(dir_plot,"DimPlot_umpa.png"),width = 10,height = 8)
  ggsave(file.path(dir_plot,"DimPlot_umpa.pdf"),width = 10,height = 8)
  
  DimPlot(in_object, reduction = "tsne",pt.size = 1)
  ggsave(file.path(dir_plot,"DimPlot_tsne.png"),width = 10,height = 8)
  ggsave(file.path(dir_plot,"DimPlot_tsne.pdf"),width = 10,height = 8)
  
  
  freq_celltype=cal_freq_NoGroup(df_label_pred,"pred_celltype_new",1)
  openxlsx::write.xlsx(freq_celltype,file.path(dir_plot,"freq_celltype.xlsx"))
  in_object
}





#For QC -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
add_qc_measurements=function(inObj,getMTgeneAuto=T,
                             MTgene=df_geneid$gene_name[grepl("^MT-",df_geneid$gene_name)],
                             ribosomalGene=df_geneid$gene_name[grepl("^RP[SL]",df_geneid$gene_name)],
                             housekeepinngGene=hkgenes_human){
  
  gene_list=rownames(x = inObj)
  
  if(getMTgeneAuto==T){
    MTgene=gene_list[grepl("^MT-",toupper(gene_list))]
    ribosomalGene=gene_list[grepl("^RP[SL]",toupper(gene_list))]
  }  
  
  # print(paste0("ribosomalGene list:",paste0(ribosomalGene,collapse = ",")))
  # print(paste0("MTgene list:",paste0(MTgene,collapse = ",")))
  # print(paste0("housekeepinngGene list:",paste0(housekeepinngGene,collapse = ",")))
  # 
  ribosomalGene1=gene_list[gene_list %in% ribosomalGene]
  MTgene1=gene_list[gene_list %in% MTgene]
  housekeepinngGene1=gene_list[gene_list %in% housekeepinngGene]
  
  print(paste0("Used ribosomalGene:",paste0(ribosomalGene1,collapse = ",")))
  print(paste0("Used MTgene:",paste0(MTgene1,collapse = ",")))
  print(paste0("Used housekeepinngGene:",paste0(housekeepinngGene1,collapse = ",")))
  
  #add qc measurement
  inObj[['percent.ribosomal']] = PercentageFeatureSet(inObj, features = ribosomalGene1)
  inObj[['percent.mt']] = PercentageFeatureSet(inObj, features = MTgene1)
  inObj[['percent.houseKeeping']] = PercentageFeatureSet(inObj, features = housekeepinngGene1)
  
  inObj
}

get_cutoff=function(inObj,
                    var_list=c("nCount_RNA","nFeature_RNA","percent.ribosomal","percent.mt","percent.houseKeeping")){
  bind_rows(lapply(var_list, function(var_in){
    x=unlist(inObj@meta.data[var_in])
    df_cutoff=data.frame(
      var=var_in,
      percentile_1=quantile(x,0.01),
      percentile_5=quantile(x,0.05),
      percentile_10=quantile(x,0.10),
      percentile_90=quantile(x,0.90),
      percentile_95=quantile(x,0.95),
      percentile_99=quantile(x,0.99),
      
      median_=median(x),
      mad_=mad(x),
      
      medianAdd3MAD=median(x)+3*mad(x),
      medianMinus3MAD=median(x)-3*mad(x),
      stringsAsFactors = F
    )
    df_cutoff
  }))
}

draw_VlnPlot_qc=function(inObj,outfile_prefix="2.data_temp/test/VlnPlot_qc"){
  
  dir_name=dirname(dirname(outfile_prefix))
  if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  dir_name=dirname(outfile_prefix)
  if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  #vilion plot
  VlnPlot1=VlnPlot(inObj, features = "nFeature_RNA", ncol = 1) + theme(legend.position = "none")
  VlnPlot2=VlnPlot(inObj, features = "nCount_RNA", ncol = 1) + theme(legend.position = "none")
  VlnPlot3=VlnPlot(inObj, features = "percent.mt", ncol = 1) + theme(legend.position = "none")
  VlnPlot4=VlnPlot(inObj, features = "percent.ribosomal", ncol = 1) + theme(legend.position = "none")
  VlnPlot5=VlnPlot(inObj, features = "percent.houseKeeping", ncol = 1) + theme(legend.position = "none")
  VlnPlot6=VlnPlot(inObj, features = "S.Score", ncol = 1) + theme(legend.position = "none")
  VlnPlot7=VlnPlot(inObj, features = "G2M.Score", ncol = 1) + theme(legend.position = "none")
  VlnPlot8=VlnPlot(inObj, features = "Singlet_Doublet_score", ncol = 1) + theme(legend.position = "none")
  
  ggarrange(VlnPlot1, VlnPlot2, VlnPlot3,VlnPlot4,VlnPlot5,VlnPlot6,VlnPlot7,VlnPlot8,ncol = 4, nrow = 2)
  ggsave(filename = paste0(outfile_prefix,".png"),width = 13,height = 9)
  ggsave(filename = paste0(outfile_prefix,".pdf"),width = 13,height = 9)
}

draw_correlationPlot_qc=function(inObj,outfile_prefix="2.data_temp/test/FeatureScatter_qc"){
  dir_name=dirname(dirname(outfile_prefix))
  if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  dir_name=dirname(outfile_prefix)
  if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  #correlation plot
  cor_plot1 <- FeatureScatter(inObj, feature1 = "nCount_RNA", feature2 = "percent.mt")+ theme(legend.position = "none")
  cor_plot2 <- FeatureScatter(inObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ theme(legend.position = "none")
  cor_plot3 = FeatureScatter(inObj, feature1 = "nCount_RNA", feature2 = "percent.ribosomal")+ theme(legend.position = "none")
  cor_plot4 = FeatureScatter(inObj, feature1 = "percent.ribosomal", feature2 = "percent.mt")+ theme(legend.position = "none")
  
  ggarrange(cor_plot1, cor_plot2, cor_plot3,cor_plot4,ncol = 4, nrow = 1)
  ggsave(filename = paste0(outfile_prefix,".png"),width = 18,height = 4)
  ggsave(filename = paste0(outfile_prefix,".pdf"),width = 18,height = 4)
}

draw_umap_qc=function(inObj,outfile_prefix="2.data_temp/test/umap_qc"){
  dir_name=dirname(dirname(outfile_prefix))
  if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  dir_name=dirname(outfile_prefix)
  if(!dir.exists(dir_name)){dir.create(dir_name)}
  
  umap1=FeaturePlot(inObj, features = "nFeature_RNA")
  umap2=FeaturePlot(inObj, features = "nCount_RNA")
  umap3=FeaturePlot(inObj, features = "percent.mt")
  umap4=FeaturePlot(inObj, features = "percent.ribosomal")
  umap5=FeaturePlot(inObj, features = "percent.houseKeeping")
  umap6=FeaturePlot(inObj, features = "S.Score")
  umap7=FeaturePlot(inObj, features = "G2M.Score")
  umap8=DimPlot(inObj,group.by='Phase')
  umap9=FeaturePlot(inObj, features = "Singlet_Doublet_score")
  umap10=DimPlot(inObj,group.by='Singlet_Doublet')
  
  ggarrange(umap1, umap2,umap3,umap4,umap5,umap6,umap7,umap8,umap9,umap10,ncol = 4, nrow = 3)
  ggsave(filename = paste0(outfile_prefix,".png"),width = 23,height = 12)
  ggsave(filename = paste0(outfile_prefix,".pdf"),width = 23,height = 12)
}


#get_cutoff_mad------------------------
get_cutoff_mad=function(inobj){
  
  Counts_per_cell = Matrix::colSums(inobj)
  Features_per_cell = Matrix::colSums(GetAssayData(inobj)>0)
  
  Counts_per_cell_upper=median(Counts_per_cell)+3*mad(Counts_per_cell)
  Counts_per_cell_lower=median(Counts_per_cell)-3*mad(Counts_per_cell)
  Counts_per_cell_lower=ifelse(Counts_per_cell_lower<0,0,Counts_per_cell_lower)
  
  Features_per_cell_upper=median(Features_per_cell)+3*mad(Features_per_cell)
  Features_per_cell_lower=median(Features_per_cell)-3*mad(Features_per_cell)
  Features_per_cell_lower=ifelse(Features_per_cell_lower<0,0,Features_per_cell_lower)
  
  out=c(Counts_per_cell_lower,Counts_per_cell_upper,Features_per_cell_lower,Features_per_cell_upper)
  names(out)=c("Counts_per_cell_lower",'Counts_per_cell_upper','Features_per_cell_lower','Features_per_cell_upper')
  out
}


Dimplot_QC=function(inObj,reduction,outfile_frefix){
  p1=FeaturePlot(object = inObj, reduction=reduction,features = "percent.mito")
  p2=FeaturePlot(object = inObj, reduction=reduction,features = "nFeature_RNA")
  p3=FeaturePlot(object = inObj, reduction=reduction,features = "nCount_RNA")
  
  p1+p2+p3
  ggsave(paste0(outfile_frefix,".png"),width = 9,height = 12)
  ggsave(paste0(outfile_frefix,".pdf"),width = 9,height = 12)
}


get_topVariableGenes=function(in_object,method_ForVariableGene="vst",N_variable){
  #Normalizing the data
  in_object <- NormalizeData(in_object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  #Identification of highly variable features (feature selection)
  in_object <- FindVariableFeatures(in_object, selection.method = method_ForVariableGene, nfeatures = N_variable)
  
  top_gene <- head(VariableFeatures(in_object), N_variable)
  
  in_object1=in_object[top_gene]
  return(in_object1)
}



# correct_batch_seurat ---------------------

batch_correction_seurat=function(inObj,var_batch,n_varaible=5000,method="RPCA"){
  cat("Avaiable methods for batch correction: RPCA, CCA, SCT\n")
  options(future.globals.maxSize= 1612185600)
  
  cat("Split object...\n")
  obj_list=SplitObject(inObj, split.by = var_batch)
  obj_list=lapply(obj_list, function(x){
    matrix_count=as.matrix(GetAssayData(x,assay = "RNA",slot = "count"))
    df_info=as.data.frame(x@meta.data)
    x_out=CreateSeuratObject(matrix_count,meta.data = df_info)
  })
  
  if(method %in% c("RPCA","CCA")){
    cat("Runing preparation\n")
    
    obj_list <- lapply(X = obj_list, FUN = function(x) {
      x <- NormalizeData(x,verbose = F)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = n_varaible,verbose = F)
    })
    
    features <- SelectIntegrationFeatures(object.list = obj_list)
    
    
    if(method=="RPCA"){
      obj_list <- lapply(X = obj_list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
      })
      cat("Finding anchors, using rPCA\n")
      obj_anchors = FindIntegrationAnchors(object.list = obj_list, anchor.features = features, reduction = "rpca",verbose = F)
    }
    
    if(method=="CCA"){
      print("Finding anchors, using CCA")
      obj_anchors=FindIntegrationAnchors(object.list = obj_list, anchor.features = features,verbose = F)
    }
    
    cat("Integrate Data\n")
    obj_combined= IntegrateData(anchorset = obj_anchors)
    
    
    DefaultAssay(obj_combined) = "integrated"
    
    obj_combined=ScaleData(obj_combined)
    obj_combined=RunPCA(obj_combined)
  }
  
  if(method =="SCT"){
    cat("Runing preparation\n")
    obj_list <- lapply(X = obj_list, FUN = function(x) {
      x <- SCTransform(x,variable.features.n = n_varaible,verbose = F)
    })
    
    features <- SelectIntegrationFeatures(object.list = obj_list,nfeatures = 10000)
    
    obj_list = PrepSCTIntegration(object.list = obj_list, anchor.features = features,verbose = FALSE)
    
    cat("Finding anchors, using rPCA\n")
    anchors_integration = FindIntegrationAnchors(object.list = obj_list, normalization.method = "SCT",
                                                 anchor.features = features, verbose = FALSE)
    
    cat("Integrate Data\n")
    obj_integrated = IntegrateData(anchorset = anchors_integration, normalization.method = "SCT",verbose = FALSE)

    obj_combined = RunPCA(obj_integrated, verbose = FALSE)
  }
  
  obj_combined
}





#SC analysis -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
SCTransform2cluster=function(objIn,resolution=1,n_pc=50,n_variableGenes=2000,n.neighbors=30){
  print("x1")
  objIn=SCTransform(objIn,variable.features.n = n_variableGenes)
  print("x2")
  objIn = RunPCA(objIn,npcs =n_pc ,verbose = FALSE)
  print("x3")
  objIn = RunUMAP(objIn, dims = 1:n_pc, verbose = FALSE,n.neighbors=n.neighbors)
  print("x4")
  objIn = FindNeighbors(objIn, dims = 1:n_pc, verbose = FALSE)
  print("x5")
  objIn = FindClusters(objIn, verbose = FALSE,resolution = resolution)
  print("x6")
  objIn
}


get_SCT_matrix=function(count_matrix,n_variableGenes=10000){
  obj_=CreateSeuratObject(count_matrix)
  
  obj_=SCTransform(obj_,variable.features.n = n_variableGenes)
  
  df_sctmatrix=as.data.frame(as.matrix(GetAssayData(obj_[["SCT"]],slot="data")))
  df_sctmatrix$gene_name=row.names(df_sctmatrix)
  df_sctmatrix1=df_sctmatrix[c("gene_name",names(df_sctmatrix)[1:(ncol(df_sctmatrix)-1)])]
  df_sctmatrix1
}


#Draw singlecell plot ---------------------------------------------------------------------------------------------------------------------------------------------------------------



























