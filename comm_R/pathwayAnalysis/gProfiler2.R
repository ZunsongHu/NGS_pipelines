#https://rdrr.io/cran/gprofiler2/f/vignettes/gprofiler2.Rmd
library(gprofiler2)

#gorth translates gene identifiers between organisms. For example, to convert gene list between mouse (source_organism = mmusculus) and human (target_organism = hsapiens):
gorth(query = c("Klf4", "Sox2", "71950"), source_organism = "mmusculus",target_organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC")

#gconvert enables to map between genes, proteins, microarray probes, common names, various database identifiers, etc, from numerous databases and for many species.
gconvert(query = c("REAC:R-HSA-3928664", "rs17396340", "NLRP1"), organism = "hsapiens",target="ENSG", mthreshold = Inf, filter_na = TRUE)

#gsnpense converts a list of SNP rs-codes (e.g. rs11734132) to chromosomal coordinates, gene names and predicted variant effects. Mapping is only available for variants that overlap with at least one protein coding Ensembl gene.
gsnpense(query = c("rs11734132", "rs4305276", "rs17396340", "rs3184504"))

gostres = gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), organism = "hsapiens")

gostplot(gostres, capped = FALSE, interactive = FALSE)

gene_list

gostres=gost(query = gene_list, organism = "mmusculus",ordered_query=T,correction_method="bonferroni")


df_enrichment=gostres$result
df_enrichment2=gostres$result


gprofiler2::(gostres)

gostres$result[1:20,]

publish_gosttable(gostres, 
                  highlight_terms = gostres$result[1:20,],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size"),
                  filename = NULL)
