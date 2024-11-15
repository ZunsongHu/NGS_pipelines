#NES > 1.5, adj.p<0.05 means significant pathway

library(clusterProfiler)
library(org.Mm.eg.db)

library(org.Hs.eg.db)

library(AnnotationDbi)


#enrichment analysis ------------------------------------------------------- ----
ego = enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)






gse_GO=clusterProfiler::enrichGO(gene  = gene_list,ont="ALL",OrgDb = "org.Mm.eg.db",keyType="SYMBOL",readable= TRUE) %>% as.data.frame()




#GSEA ------------------------------------------------------- ----
df_DE_2=df_DE_1 %>% arrange(desc(avg_log2FC))
geneList=df_DE_2$avg_log2FC
names(geneList)=df_DE_2$gene
# geneList[1:10]
cat("Running gseGO...\n")
out_gseGO = gseGO(geneList= geneList,OrgDb= "org.Mm.eg.db",ont= "ALL",keyType="SYMBOL",pvalueCutoff = 0.05)
df_gseGO=out_gseGO %>% as.data.frame()
write_tsv(df_gseGO,paste0(prefix,".df_GSEgo",".tsv"))































