library(dplyr)
library(GEOquery)
gse="gse54618"


get_gseData=function(gse,get_exp=T){
  
  if(get_exp){
    gset = getGEO(gse, GSEMatrix =T, getGPL=T)
    gset1=gset[[1]]
    
    df_coldata=pData(gset1) %>% as.data.frame()
    df_anno=fData(gset1)
    df_exp=exprs(gset1) %>% as.data.frame()
    return(list(coldata=df_coldata,exp=df_exp,anno=df_anno))
  }
  
  if(!get_exp){
    gset = getGEO(gse, GSEMatrix =F, getGPL=F)
    gset1=gset[[1]]
    
    df_coldata=pData(gset1) %>% as.data.frame()
    return(list(coldata=df_coldata))
  }
  gc()
}

get_geneMean_fromProbeID=function(exp,gene.anno,var_id,var_gene){
  gene.anno1=gene.anno[c(var_id,var_gene)]
  names(gene.anno1)=c("ID_REF","gene")
  #
  gene.anno1<-gene.anno1[!is.na(gene.anno1$gene),]
  
  
  if(!'ID_REF' %in% colnames(exp)){exp$ID_REF=row.names(exp)}
  # exp$ID_REF=row.names(exp)
  # if(!all(exp$ID_REF %in% gene.anno1$ID_REF)){stop("row name of exp not in anno names")}
  
  gse.dat<-merge(gene.anno1,exp,by='ID_REF')
  gse.dat=gse.dat[,-1]
  
  #calculation mean
  gse.dat=aggregate(.~gene,gse.dat,mean)
  rownames(gse.dat)=gse.dat$gene
  gse.dat1=gse.dat[,-1]
  
  gse.dat1
}


# extract exp and clinical data -------------------------------------------------- ----
setwd("c:/Users/zuhu/OneDrive - City of Hope National Medical Center/project/211/2024.10.15_eclampsia/")
gse_id="GSE54618"
file_out="0.original/GSE54618.rds"

# if(!dir.exists(dirname(file_out))){dir.create(dirname(file_out))}


obj_gse=get_gseData(gse_id)
# obj_gse$exp
obj_gse$exp[1:10,1:10]
df_anno_=obj_gse$anno

obj_gse$anno1=obj_gse$anno %>%
  mutate(gene_symbol=Symbol) %>%
  dplyr::select(ID,gene_symbol)

obj_gse$exp_gene=get_geneMean_fromProbeID(exp = obj_gse$exp,gene.anno = obj_gse$anno1,var_id = "ID",var_gene = "gene_symbol")

df_col_=obj_gse$coldata

table(df_col_$source_name_ch1)

obj_gse$df_col=obj_gse$coldata %>%
  transmute(id=geo_accession,
            group=case_when(
              source_name_ch1=="normotensive pregnancy" ~ "normal",
              source_name_ch1=="pregnancy complicated with preeclampsia" ~ "preeclampsia",
              source_name_ch1=="pregnancy complicated with preeclampsia and HELLP syndrome" ~ "preeclampsia",
              TRUE~NA
            ),
            source_name_ch1=source_name_ch1)

table(obj_gse$df_col$group)

saveRDS(obj_gse,file_out)


# group=case_when(
#   characteristics_ch1.1=="disease state: control" ~ "control",
#   characteristics_ch1.1=="disease state: PMF patient" ~ "PMF",
#   TRUE~NA
# )









