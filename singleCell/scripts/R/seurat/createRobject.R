

singleR1=readRDS("analysis/obj/A11419.mouse_9_main.SingleR_out.rds")
singleR2=readRDS("analysis/obj/A11420.mouse_9_main.SingleR_out.rds")

obj_1=readRDS("analysis/obj/A11419.SeuratObj.QC.rds")
obj_2=readRDS("analysis/obj/A11420.SeuratObj.QC.rds")


get_obj_new=function(obj_,singleR_){
  df_matrix=GetAssayData(obj_,assay = "RNA",layer = "counts")
  df_meta=obj_@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("barcode") %>% 
    left_join(singleR_$df_singleR_out %>% mutate( seurat_clusters=NULL))
  row.names(df_meta)=df_meta$barcode
  
  obj_new=CreateSeuratObject(counts = df_matrix,meta.data = df_meta)
  obj_new
}
