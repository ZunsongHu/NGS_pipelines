Idents(obj_)="cd45_g"


df_DE_=FindMarkers(obj_,ident.1 = "CD45.2",ident.2 = "CD45.1",logfc.threshold=0,min.pct=0)
df_DE_1=df_DE_ %>% tibble::rowid_to_column(var="gene")

write_tsv(df_DE_1,paste0(prefix,".df_DE",".tsv"))

