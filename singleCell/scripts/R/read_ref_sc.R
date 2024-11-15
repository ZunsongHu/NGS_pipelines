dir_pipeline="/net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/"

# "/net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/ref/grch38/gtf/V102/Homo_sapiens.GRCh38.V102.summary.txt"
# "/net/nfs-irwrsrchnas01/labs/zgu_grp/Individual/zuhu/pipeline/ref/cellMarkers/HouseKeepingGenes_mouse.txt"

dir_ref=paste0(dir_pipeline,"ref/")

#get ENSG id ---------------------------------------------------------------------------------------------------------------------------------------------------------------

hkgenes_human=readLines(paste0(dir_ref,"cellMarkers/HouseKeepingGenes.txt"))
hkgenes_mouse=readLines(paste0(dir_ref,'cellMarkers/HouseKeepingGenes_mouse.txt'))

file_all_geneid=paste0(dir_ref,"grch38/gtf/V102/Homo_sapiens.GRCh38.V102.summary.txt")
# geneid_all=read.table(file_all_geneid,header = T,stringsAsFactors = F,sep="\t")

df_geneid=read.table(file_all_geneid,header = T,stringsAsFactors = F)

get_ENSGID=function(genes){x1=df_geneid[df_geneid$gene_name %in% genes,];x1$gene_id}

s_ENSG=get_ENSGID(Seurat::cc.genes$s.genes)
g2m_ENSG=get_ENSGID(Seurat::cc.genes$g2m.genes)
MT_ENSG=get_ENSGID(df_geneid$gene_name[grepl("^MT-",df_geneid$gene_name)])
Ribosomal_ENSG=get_ENSGID(df_geneid$gene_name[grepl("^RP[SL]",df_geneid$gene_name)])
housekeeping_ENSG=get_ENSGID(hkgenes_human)

rm(file_all_geneid)
rm(dir_pipeline)
rm(dir_ref)





















