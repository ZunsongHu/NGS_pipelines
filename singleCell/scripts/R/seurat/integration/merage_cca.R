
### 使用CAA方法整合数据（其它方法见前期教程）
### 对Seurat对象进行切割
seu_list <- SplitObject(seu_merge,split.by="orig.ident")
### 进行均一化处理寻找高变基因
for (i in 1:length(seu_list)) {  
  seu_list[[i]] <- NormalizeData(seu_list[[i]])  
  seu_list[[i]] <- FindVariableFeatures(seu_list[[i]], selection.method = "vst",nfeatures = 3000)}

### 使用多个核进行单细胞数据分析library(future)plan("multiprocess", workers =3) 

###计划用几个线程跑
options(future.globals.maxSize = 2000 * 1024^2)
###给每个线程分配多少G的运行内存1024^2为1M,表示为一线程2

Gfeatures <- SelectIntegrationFeatures(object.list = seu_list) 
###选择整合基因也可不选择###寻找锚定点

seu.anchors <- FindIntegrationAnchors(object.list = seu_list,anchor.features = features)

###根据锚定点整合数据

seu_merge <- IntegrateData(anchorset = seu.anchors)

### 将activate assay更改为integrated，integrated用于后续FindClusters，和作图

DefaultAssay(seu_merge) <- "integrated"
seu_merge=ScaleData(seu_merge)
seu_merge <- RunPCA(seu_merge, npcs = 30, verbose = T)
seu_merge <- FindNeighbors(seu_merge, reduction = "pca", dims = 1:20)
seu_merge <- FindClusters(seu_merge, resolution = 0.8)

### 进行Umap和Tsne的可视化

seu_merge <- RunUMAP(seu_merge, reduction = "pca", dims = 1:20)
plot1 = DimPlot(seu_merge, reduction = "umap", group.by='orig.ident',label.size = 0.5,pt.size = 0.1,cols = colors,raster.dpi=c(1024,1024))+  
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        legend.position = "right")+  
  labs(title = "Sample Umap")

plot1