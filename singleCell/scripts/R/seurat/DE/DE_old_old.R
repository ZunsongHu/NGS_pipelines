
Idents(obj_subset)="group"

obj_subset=PrepSCTFindMarkers(obj_subset)

df_DE_macrophage=FindMarkers(obj_subset,ident.1 = "ko",ident.2 = "wt")

df_DE_macrophage1=df_DE_macrophage %>% tibble::rownames_to_column(var = 'feature') %>%  mutate(contrast="KO_vs_wt")
write_tsv(df_DE_macrophage1,"analysis/DE/ko_vs_wt.in_antiTumorMacrophage/df_DE_antiTumormacrophage.tsv")

gg_volcanoPlot(df_in = df_DE_macrophage1,P = "p_val_adj",FC = "avg_log2FC",var_label = "feature",n_top_to_show = 20,prefix = "analysis/DE/ko_vs_wt.in_antiTumorMacrophage/volcanoPlot",w = 8,h = 6)

















