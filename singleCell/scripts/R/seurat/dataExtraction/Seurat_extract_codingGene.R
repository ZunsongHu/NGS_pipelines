# rm(list=ls())
source("/home/zgu_labs/bin/R/SingleCell/F_sc.R")

# #temp parameters ----
# setwd("/scratch/zuhu/project/srivi/singlecell/")
# rds_in="out_raw/health_5p5K/Seurat/health_5p5K.SeuratObj.raw.rds"
# rds_out="test/test.rds"
# species="mouse"

#get parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

rds_in=args[1]
rds_out=args[2]
species=args[3]

#create dir out ----
if(!dir.exists(dirname(rds_out))){dir.create(dirname(rds_out))}

#read rds
cat("reading obj ...\n")
obj_raw=readRDS(rds_in)

cat("Raw object:\n")
obj_raw
cat("\n")

names_=row.names(x=obj_raw)

#get coding gene names ----
if(species=="mouse"){
  info_gtf=read_tsv("/ref_genomes/genome_anno/mouse/gtf/V102/Mus_musculus.GRCm38.V102.withChr.summary.txt")
}

if(species=="human"){
  info_gtf=read_tsv("/ref_genomes/genome_anno/human/gtf/V102/Homo_sapiens.GRCh38.V102.summary.txt")
}

coding_genes=info_gtf$gene_name[info_gtf$gene_biotype=="protein_coding"]

#get overlap ----

overlap_=names_[names_ %in% coding_genes]

obj_raw1=obj_raw[overlap_,]

cat("New object:\n")
obj_raw1
cat("\n")
cat("Saving output object ...\n")
saveRDS(obj_raw1,rds_out)


cat("Done\n")

