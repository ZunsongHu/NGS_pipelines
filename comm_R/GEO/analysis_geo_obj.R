

gc(rm(list = ls()))

file_geo_list="0.original/list_gse24129.rds"
dir_out="analysis/gse24129/"
id="gse24129"
comparisons=c("pre_eclampsia_vs_normal","FGR_vs_normal")

file_geo_list="0.original/GSE54618.rds"
dir_out="analysis/gse54618/"
id="gse54618"
comparisons=c("preeclampsia_vs_normal")

file_geo_list="0.original/GSE74341.rds"
dir_out="analysis/GSE74341/"
id="GSE74341"
comparisons=c("EOPE_vs_atTerm","LPOE_vs_atTerm","EOPE_vs_LPOE")

dir.create("")




log2_transfrom=F
scale_=F
w_boxplot=10
h_boxplot=5

w_pheatmap=12
h_pheatmap=8

w_pca=12
h_pca=8

w_gsea=10
h_gsea=10

cutoff_quantile_exp=0.5

# p_DE="P.Value"
p_DE="adj.P.Val"
n_top_to_show_volcano=10
P_cutoff=0.05
FC_cutoff=0.5

n_top_DE=20
featuresToShow="DIAPH3"

# library ------------------------------------------------------------ ----
library(pheatmap)
library(dplyr)
library(ggplot2)
library(limma)
library(stringr)

library(clusterProfiler)
library(org.Hs.eg.db)

library(ggalign)

db_ref="org.Hs.eg.db"



# functions ---------------------------------------------------------------- ----
run_DE_limma=function(exp,df_info,var_name,comparisons,covariates=NULL,adjust_method="BH"){
  if (!var_name %in% colnames(df_info)) {stop(paste("Variable", var_name, "not found in df_info"))}
  
  comparison_levels=c(word(comparisons,1,sep="_vs_"),word(comparisons,2,sep="_vs_"))
  if(!all(comparison_levels %in% df_info[[var_name]])){stop("comparison level not in data")}
  
  #get df_info for analysis
  df_info_subset <- df_info[df_info[[var_name]] %in% comparison_levels, ]
  df_info_subset[[var_name]] <- factor(df_info_subset[[var_name]], levels = comparison_levels)
  
  #get exp for analysis 
  common_samples = intersect(colnames(exp), rownames(df_info_subset))
  if (length(common_samples) == 0) {stop("No matching samples between expression data and df_info")}
  exp_subset <- exp[, common_samples, drop = FALSE]
  df_info_subset <- df_info_subset[common_samples, , drop = FALSE]
  
  #get formula
  formula_terms <- c("0", var_name)
  if (!is.null(covariates)) {formula_terms <- c(formula_terms, covariates)}
  
  design_formula <- as.formula(paste("~", paste(formula_terms, collapse = " + ")))
  design = model.matrix(design_formula, data = df_info_subset)
  colnames(design) <- sort(levels(df_info_subset[[var_name]]))
  
  #run lmFit
  fit <- lmFit(exp_subset, design)
  
  #get contrast matrix
  contrast_formula <- paste0(comparison_levels[1], " - ", comparison_levels[2])
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  
  #run contrast
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  #get DE table
  results <- topTable(fit2, number = Inf, adjust.method = adjust_method) %>% 
    tibble::rownames_to_column(var="gene") %>% 
    mutate(contrast=comparisons)
  results
}


# QC ---------------------------------------------------------------- ----
{
  # read obj ----------------------------------------------------------------- ----
  obj_list=readRDS(file_geo_list)
  
  if(!dir.exists(dir_out)){dir.create(dir_out)}
  
  # get Exp  ------------------------------------------------------------ ----
  exp=obj_list$exp_gene
  df_col=obj_list$df_col
  
  exp=exp[,df_col$id]
  
  # Log Transformation  ------------------------------------------------------------ ----
  if(log2_transfrom){exp=log2(exp+1)}
  
  # standardizaton (Z-score normalization)  ------------------------------------------------------------ ----
  #Centers the data around zero with a standard deviation of one.
  if(scale_){exp=scale(exp)}
  
  # summary  ------------------------------------------------------------ ----
  sink(paste0(dir_out,id,".log"))
  summary(exp)
  sink()
  
  # boxPlot for exploration  ------------------------------------------------------------ ----
  pdf(paste0(dir_out,"boxplot.exp.",id,".pdf"),width = w_boxplot,height = h_boxplot)
  boxplot(exp,outline=FALSE)
  dev.off()
}


# Run analysis ---------------------------------------------------------------- ----
{
  # correlation ------------------------------------------------------------ ----
  corMatrix <- cor(exp,use="c")
  
  sample_info=obj_list$df_col["group"]
  row.names(sample_info)=obj_list$df_col$id
  
  pdf(paste0(dir_out,"pheatmap.exp.",id,".pdf"),width = w_pheatmap,height = h_pheatmap)
  pheatmap(corMatrix,annotation_col=sample_info)   
  dev.off()
  
  # pca ------------------------------------------------------------ ----
  pca <- prcomp(t(exp))
  
  all(row.names(sample_info)==row.names(pca$x))
  
  cbind(sample_info, pca$x) %>% 
    ggplot(aes(x = PC1, y=PC2, col=group)) + geom_point() 
  ggsave(paste0(dir_out,"pcaPlot.",id,".pdf"),width = w_pca,height = h_pca)
  
  # DE analysis  ------------------------------------------------------------ ----
  table(obj_list$df_col$group)
  
  # compare=comparisons[1]
  
  for(compare in comparisons){
    #run DE
    if(!file.exists(paste0(dir_out,"df_DE.",compare,".xlsx"))){
      df_DE=run_DE_limma(exp = exp,df_info = df_col,var_name = "group",comparisons = compare)
      df_DE %>% openxlsx::write.xlsx(paste0(dir_out,"df_DE.",compare,".xlsx"))
    }
    if(file.exists(paste0(dir_out,"df_DE.",compare,".xlsx"))){
      df_DE=readxl::read_excel(paste0(dir_out,"df_DE.",compare,".xlsx"))
    }
    
    #volcano plot
    gg_volcanoPlot(df_in = df_DE,P = p_DE,FC = "logFC",var_label = "gene",P_cutoff = P_cutoff,FC_cutoff = FC_cutoff,n_top_to_show = n_top_to_show_volcano,
                   featuresToShow = featuresToShow,show_DE_number = T,
                   prefix = paste0(dir_out,"df_DE.",compare),w = 6,h = 4,write_df_DE = F)
    
    #DE heatmap
    df_DE1=df_DE %>% mutate(g=ifelse(logFC>0,"Pos",ifelse(logFC<0,"Neg",NA))) %>% 
      group_by(g) %>% arrange(P.Value) %>% slice_head(n=n_top_DE)
    
    df_anno_col=df_col["group"] %>% tibble::rownames_to_column(var="id") %>% 
      filter(group %in% c(word(compare,1,sep="_vs_"),word(compare,2,sep="_vs_")))
    df_anno_col=df_anno_col[order(df_anno_col$group),]
    exp1=exp[df_DE1$gene,df_anno_col$id]
    row.names(df_anno_col)=df_anno_col$id
    df_anno_col$id=NULL
    
    pdf(paste0(dir_out,"df_DE.",compare,".pheatmap.pdf"),w=10,h=10)
    pheatmap::pheatmap(exp1,annotation_col = df_anno_col,cluster_cols = F)
    dev.off()
    
    #GSEA
    df_DE1=df_DE %>% arrange(desc(logFC))
    
    geneList=df_DE1$logFC
    names(geneList)=df_DE1$gene
    
    if(!file.exists(paste0(dir_out,"out_gseGO.",compare,".rds"))){
      out_gseGO = gseGO(geneList= geneList,OrgDb= db_ref,ont= "ALL",keyType="SYMBOL",pvalueCutoff = 0.05)
      saveRDS(out_gseGO,paste0(dir_out,"out_gseGO.",compare,".rds"))
    }
    if(file.exists(paste0(dir_out,"out_gseGO.",compare,".rds"))){
      out_gseGO=readRDS(paste0(dir_out,"out_gseGO.",compare,".rds"))
    }
    
    df_gseGO=out_gseGO %>% as.data.frame()
    df_gseGO %>% openxlsx::write.xlsx(paste0(dir_out,"df_gseGO.",compare,".xlsx"))
    
    gg_gseaNES_clusterProfiler(df_gsea = df_gseGO,topN = 20,prefix = paste0(dir_out,"out_gseGO.",compare),w = w_gsea,h = h_gsea,title = compare)
    
  }
}

