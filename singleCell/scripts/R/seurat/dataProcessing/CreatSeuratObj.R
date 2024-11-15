source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

#temp parameters
# setwd("/scratch/zuhu/project/xiaojingYan/")
# # 
# dir_matrix="out_raw/test1/Cellranger_Count_Expectnumber/filtered_feature_bc_matrix/"
# id="test1"
# outfile="analysis/obj/test1.SeuratObj.raw.rds"

#get parameters
args <- commandArgs(trailingOnly = TRUE)
print(args)

dir_matrix=args[1]
id=args[2]
outfile=args[3]

#create outfile folder
if(!dir.exists(dirname(outfile))){dir.create(dirname(outfile))}

#create seurat object 
seurat_obj=create_seurat_object(dir_matrix,
                     project_name=id,
                     min.cells=3,
                     min.features=0,
                     gene.column=2)

seurat_obj[["group"]]=id
print(seurat_obj)

saveRDS(seurat_obj,file=outfile)

# 
# df_barcode=get_barcode_df(seurat_obj) %>% sample_n(300)
# 
# seurat_obj_test=subset(seurat_obj,cells=df_barcode$barcode)
# 
# saveRDS(seurat_obj_test,file=
#           "/scratch/zuhu/project/ZhaohuiGu/singlecell/PublicData/2020.GenomeMedicine.ETV6RUNX1/out/scRNA_BMMC_D1T1_Expression1/Seurat/test.rds")
# 



