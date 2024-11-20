obj_list=list()
obj_list[[i]]=obj_new
obj_all <- merge(obj_list[[1]], 
                 y = c(obj_list[[2]],obj_list[[3]],obj_list[[4]],obj_list[[5]],obj_list[[6]],obj_list[[7]],obj_list[[8]],obj_list[[9]],obj_list[[10]],obj_list[[11]]),
                 add.cell.ids = c("A2","A3","A4","A5","C1","C3","C4","C5","health_10K","health_5K","health_v2"), 
                 project = "meged")



#merge 2
obj_all <- merge(obj1_34.1, 
                 y = c(obj1_34.2,obj1_14.1,obj1_14.2),
                 add.cell.ids = c("CD34.1","CD34.2","CD14.1","CD24.2"), 
                 project = "meged")



#log normalization
obj_all <- NormalizeData(obj_all)
obj_all <- FindVariableFeatures(obj_all)
obj_all <- ScaleData(obj_all)
obj_all <- RunPCA(obj_all)
obj_all <- RunUMAP(obj_all, dims = 1:50)

#cca merge
obj_all <- IntegrateLayers(object = obj_all, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
    verbose = FALSE)
obj_all[["RNA"]] <- JoinLayers(obj_all[["RNA"]])
obj_all <- FindNeighbors(obj_all, reduction = "integrated.cca", dims = 1:50)
obj_all <- FindClusters(obj_all, resolution = 1)
obj_all <- RunUMAP(obj_all, dims = 1:50, reduction = "integrated.cca")

#SCT
obj_all <- SCTransform(obj_all)
obj_all <- RunPCA(obj_all)
obj_all <- RunUMAP(obj_all, dims = 1:50)

#SCT merge
obj_all <- IntegrateLayers(object = obj_all, method = CCAIntegration, normalization.method = "SCT", verbose = F)
obj_all <- RunUMAP(obj_all, dims = 1:50, reduction = "integrated.dr")
obj_all <- FindNeighbors(obj_all, reduction = "integrated.dr", dims = 1:50)
obj_all <- FindClusters(obj_all, resolution = 1)


#for two samples -----------------------------------------------------------------------------------------------------------
obj_wt=get_obj_new(obj_1,singleR1)
obj_ko=get_obj_new(obj_ = obj_2,singleR_ = singleR2)
obj_wt$group="wt"
obj_ko$group="ko"
obj_merge <- merge(obj_wt,
                   y = c(obj_ko),
                   add.cell.ids = c("wt","ko"),
                   project = "yuanChe")

obj_merge = SCTransform(obj_merge)
obj_merge = RunPCA(obj_merge)
obj_merge = IntegrateLayers(object = obj_merge, method = CCAIntegration, normalization.method = "SCT", verbose = F)
obj_merge = RunUMAP(obj_merge, dims = 1:50, reduction = "integrated.dr")
obj_merge = FindNeighbors(obj_merge, reduction = "integrated.dr", dims = 1:50)
obj_merge = FindClusters(obj_merge, resolution = 1)
saveRDS(obj_merge,"analysis/obj/obj_merge.rds")




