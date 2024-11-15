# source("/net/nfs-irwrsrchnas01/labs/zgu_grp/DryLab/bin/R_apollo/SingleCell/F_sc.R")

#usage: 
# Rscript_seurat5 /net/nfs-irwrsrchnas01/labs/zgu_grp/DryLab/bin/R_apollo/SingleCell/gsea/run_gsea_escape.R \
# analysis/obj/obj.11samples.logNormalizationBatchCorrected.selected.rds analysis/gsea/ \
# Bax_Bak1_Bok_Bad_Bid_Noxa1_Bik_Bmf_Hrk_Bnip3_Bcl2l11_Bbc3 \
# proApoptosis_BCL2family

#load libraries ---------------------------------------------------------------------------------------------------------------------
options(warn=-1)
suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(escape)
})
options(warn=0)

#temp parameters  ---------------------------------------------------------------------------------------------------------------------
# file_obj="/net/nfs-irwrsrchnas01/labs/zgu_grp/Group/Grp_Swaminathan/singleCell/analysis/obj/obj.11samples.logNormalizationBatchCorrected.selected.rds"
# dir_out="/net/nfs-irwrsrchnas01/labs/zgu_grp/Group/Grp_Swaminathan/singleCell/analysis/gsea/"
# genes_in_pathway="Bax_Bak1_Bok_Bad_Bid_Noxa1_Bik_Bmf_Hrk_Bnip3_Bcl2l11_Bbc3"
# name_pathway="proApoptosis_BCL2family"

#arguments  ---------------------------------------------------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
print(args)
file_obj=args[1]
dir_out=args[2]
genes_in_pathway=args[3]
name_pathway=args[4]

#phase parameters ---------------------------------------------------------------------------------------------------------------------
if(!dir.exists(dir_out)){dir.create(dir_out)}

gene.sets <- list(pathway = c(unlist(strsplit(genes_in_pathway,split = "_"))))
names(gene.sets)=name_pathway

file_out=paste0(dir_out,"df_score_escape.",name_pathway,'.tsv')

#functions  ---------------------------------------------------------------------------------------------------------------------
write_tsv=function(indata,outfile){
  if(!dir.exists(dirname(outfile))){dir.create(dirname(outfile))}
  write.table(indata,outfile,col.names = T,row.names = F,quote = F,na = "",sep="\t")
}

read_tsv=function(infile){
  read.table(infile,header = T,sep = "\t",stringsAsFactors = F)
}

#read obj  ---------------------------------------------------------------------------------------------------------------------
cat("Reading obj...\n")

obj_=readRDS(file_obj)

# set.seed(100)
# sampled_cells =sample(Cells(obj_),5000)
# obj_ = subset(obj_, cells = sampled_cells)

# calculate the enrichment score   ---------------------------------------------------------------------------------------------------------------------
# method_="GSVA"
# method_="ssGSEA"
# 
# method_="AUCell"
# method_="UCell"
cat("Calculate escape score...\n")
get_escape_values=function(obj_,method_="AUCell"){
  obj_ = runEscape(obj_, 
                    method = method_,
                    gene.sets = gene.sets, 
                    groups =  2000, 
                    min.size = 5,
                    make.positive=T,
                    new.assay.name = "escape")
  
  escape_values=GetAssayData(obj_, assay = "escape")[1:ncol(obj_)]
  escape_values
}

df_escape=data.frame(barcode=row.names(obj_@meta.data),stringsAsFactors = F) %>% 
  mutate(
    score_escape_AUCell=get_escape_values(obj_,method_="AUCell"),
    score_escape_UCell=get_escape_values(obj_,method_="UCell"),
  )

write_tsv(df_escape,file_out)

cat(paste0("Done escape for ",name_pathway,"!\n"))

























