library(VennDiagram)

file_="analysis/gse24129/df_DE.FGR_vs_normal.xlsx"

get_featureList=function(file_){
  df1=readxl::read_excel(file_) %>% 
    filter(logFC>0) %>% arrange(P.Value) %>% 
    slice_head(n=200)
  
  df1$gene
}


set1=get_featureList("analysis/gse24129/df_DE.FGR_vs_normal.xlsx")
set2=get_featureList("analysis/gse24129/df_DE.pre_eclampsia_vs_normal.xlsx")
set3=get_featureList("analysis/gse54618/df_DE.preeclampsia_vs_normal.xlsx")
set4=get_featureList("analysis/GSE74341/df_DE.EOPE_vs_atTerm.xlsx")
set5=get_featureList("analysis/GSE74341/df_DE.EOPE_vs_LPOE.xlsx")

cols_4=c('red4','skyblue','#CCCC33','magenta3',)
cols_5=c('#DE6E67','#717335',"#2AB47A",'#34ACE2',"#AF70A5")

p1=venn.diagram(
  x = list(Set1 = set1, Set2 = set2, Set3 = set3,Set4=set4,Set5=set5),
  category.names = c("up in FGR", "up in pre_eclampsia(gse24129)", "up in preeclampsia(gse54618)","up in EOPE_vs_atTerm","up in EOPE_vs_LPOE"),
  filename = NULL,  
  output = TRUE,
  fill = cols_5,
  alpha = 0.5,
  cat.pos = 0,
  cat.dist = 0.05,
  margin = 0.1
)

pdf("analysis/venn.diagram.pdf",width=8,height = 8)
grid.draw(p1)
dev.off()


