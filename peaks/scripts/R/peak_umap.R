# set ---------------------------
source("/home/zgu_labs/bin/R/functions.R")
source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

# testing parameters ----
setwd("/scratch/zuhu/project/Grp_NadiaCarlesso/atac/")

# dir_out=        "test/"
# label=          "test"
# file_obj=       "analysis/second_23samples/obj/obj.second_23samples.peak.cpm.rds"
# var_group=      "group"
# feature_panel=  "MAD1000"
# variable_n=     500
# assay_name_in=  "cpm"

# parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

dir_out=        args[1]
label=          args[2]
file_obj=       args[3]
var_group=      args[4]
feature_panel=  args[5]
variable_n=     as.numeric(args[6])
assay_name_in=  args[7]

# create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

# read obj ----
message("read")
obj_=readRDS(file_obj)

df_info=as.data.frame(colData(obj_$SE))

#run umap ----
message("umap")
obj_=run_umap(obj_in = obj_,n_neighbors = 10,variable_n = variable_n,feature_panel = feature_panel,assay_name_in = assay_name_in)

#draw umap ----
message("draw umap")
group_levels=unique(unlist(df_info[var_group]))

col_in=cols_in_default[1:length(group_levels)]
names(col_in)=group_levels

draw_DimPlot(obj_,group.by = var_group,cols = col_in,reduction = "umap",size = 2)
ggsave(paste0(dir_out,"umap.",label,".variableN",variable_n,".png"),height = 4,width = 5)
ggsave(paste0(dir_out,"umap.",label,".variableN",variable_n,".pdf"),height = 4,width = 5)

#save umap ----
message("save output")
df_umap=get_embeding_feature(obj_in = obj_,features = c(var_group),assay_name_in = "cpm",reduction = "umap")
write_tsv(df_umap,paste0(dir_out,"df_umap.",label,".variableN",variable_n,".tsv"))























