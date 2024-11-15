# Gene list in --------------------------------
# file.copy("//smb-irwrsrchnas01/labs/zgu_grp/DryLab/refData/genome_anno/human/gtf/V102/Homo_sapiens.GRCh38.V102.summary.txt","ref/Homo_sapiens.GRCh38.V102.summary.txt")
# file.copy("//smb-irwrsrchnas01/labs/zgu_grp/DryLab/bin/R/RNAseq/GEP_prediction/DUX4_gene.list","ref/DUX4_gene.list")
# file.copy("//smb-irwrsrchnas01/labs/zgu_grp/DryLab/refData/genome_anno/mouse/gtf/V102/Mus_musculus.GRCm38.V102.withChr.summary.txt","ref/Mus_musculus.GRCm38.V102.withChr.summary.txt")

dir_ref="c:/Users/zuhu/OneDrive - City of Hope National Medical Center/project/ref/"

file_gene_infor=paste0(dir_ref,"Homo_sapiens.GRCh38.V102.summary.txt")
dux4_genelist=paste0(dir_ref,"DUX4_gene.list")

info_gtf=read.table(file_gene_infor,header = T,stringsAsFactors = F)
codinggeneList=info_gtf$gene_id[info_gtf$gene_biotype=="protein_coding" | info_gtf$gene_id %in% dux4_genelist]
geneXYM=info_gtf$gene_id[info_gtf$chr %in% c('chrM','chrX','chrY')]
gene_in=codinggeneList[!codinggeneList %in% geneXYM]

info_gtf_hg38=info_gtf

info_gtf1=info_gtf %>% filter(gene_id %in% codinggeneList) %>% 
  group_by(gene_name) %>% dplyr::mutate(n=n()) %>% arrange(desc(n),gene_id) %>% 
  slice_head(n=1)

codinggeneList_noDup=info_gtf1$gene_id

gtf_mouse=read_tsv(paste0(dir_ref,"Mus_musculus.GRCm38.V102.withChr.summary.txt"))