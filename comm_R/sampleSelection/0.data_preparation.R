
gc(rm(list=ls()))

source("/home/zgu_labs/bin/R/RNAseq/F_RNA.R")

#parameters ----
# args <- commandArgs(trailingOnly = TRUE)
# print(args)
# 
# vst_rds=      args[1]
# dataset_label=args[2]
# file_batch=   args[3]
# var_batch=    args[4]
# dir_out=      args[5]

#parameters ----
setwd("/scratch/zuhu/project/BALL/tsne/")

#read id ----
df_info_2955=readxl::read_excel("id/Master_BALL_2955Samples.xlsx")

paste0(names(df_info_2955),collapse = "','")

table(df_info_2955$library_1)

df_info_2955_1=df_info_2955[c('COH_sample','sample_id','Disease','library','fusion','PrimarySubtype','cohort','ReadNumber','MapRate','DupRates','Depth','Coverge_ge_30X',
                            'PAX5.mut','IL7R.mut','PTPN11.mut','NF1.mut','TP53.mut','NRAS.mut','ZEB2.mut','IKZF1.mut','JAK3.mut','FLT3.mut','KRAS.mut','DNMT3A.mut',
                            'European','Native','Eastasian','Southasian','African',
                            'B_percentage',
                            'sex_RNAseqCNV','CHR.N_RNAseqCNV',
                            'ZEBP2','NUTM1','IKZF1','PAX5','CRLF2','DUX4','DUX4L6',
                            'batch1','batch2','batch3','batch3_1','library_1')] %>%
  mutate(diag=PrimarySubtype)

write_tsv(df_info_2955_1,"id/df_info_2955_1.tsv")

table(df_info_2955_1$diag)

draw_diag_barplot(df_info_2955_1,"sample_selection/2955/diag_barplot_2955.pdf")

df_for_prediction=df_info_2955_1 %>%
  filter(!diag %in% c("CRLF2(non-Ph-like)","ETV6-RUNX1-like","KMT2A-like","Ph-like","ZNF384-like","Near haploid","Low hyperdiploid","Other"))

write_tsv(df_for_prediction,"id/df_info_2391.tsv")

#read vst ----
vst_2239=read_tsv("vst_raw/vst_gene_2239.tsv")
vst_TargetALL=read_tsv("vst_raw/TargetALL.tsv")
vst_PanLeukemia1=read_tsv("vst_raw/Pan_Leukemia566.tsv")
vst_PanLeukemia2=read_tsv("vst_raw/Pan_Leukemia227.tsv")
load("vst_raw/vst_gene_EGA_432samples.rdata")

vst_all=vst_2239 %>% 
  left_join(vst_TargetALL) %>% 
  left_join(vst_PanLeukemia1) %>% 
  left_join(vst_PanLeukemia2) %>%
  left_join(df_vst_EGA)

names(vst_all)[1:10]

saveRDS(vst_all,"vst_raw/vst_gene_4350.rds")

vst_2955=vst_all[c("feature",df_info_2955$COH_sample)]
row.names(vst_2955)=vst_2955$feature

saveRDS(vst_2955,"vst_raw/vst_gene_2955.rds")

vst_2955=readRDS("vst_raw/vst_gene_2955.rds")

vst_2391=vst_2955[c("feature",df_for_prediction$COH_sample)]

saveRDS(vst_2391,"vst_raw/vst_gene_2391.rds")










