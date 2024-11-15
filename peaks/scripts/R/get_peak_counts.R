# set ---------------------------
library(GenomicRanges)
library(csaw)
source("/home/zgu_labs/bin/R/functions.R")
# testing parameters ----
# setwd("/scratch/zuhu/project/Grp_NadiaCarlesso/atac/")
# 
# dir_out=        "test/"
# label=          "test"
# file_info=      "analysis/obj/df_info_24samples.tsv"
# file_peak=      "peak.23.list"
# file_bam=       "bam.23.list"
# file_blacklist= "/ref_genomes/genome_anno/mouse/mm10-blacklist.v2.bed"
# species=        "mouse"

# parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

dir_out=        args[1]
label=          args[2]
file_info=      args[3]
file_peak=      args[4]
file_bam=       args[5]
file_blacklist= args[6]
species=        args[7]

# create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

# read peaks ----
files_peak=readLines(file_peak)

message("\n")
df_peaks=bind_rows(lapply(files_peak, function(file_peak_one){
  df_one=read.table(file_peak_one,sep="\t",stringsAsFactors = F)[,1:3] %>% mutate(file=file_peak_one)
  names(df_one)[1:3]=c("chrom","start","end")
  message(paste0("  ",nrow(df_one)," peaks found in ",file_peak_one,"."))
  df_one
})) 

message(paste0("\nTotally ",nrow(df_peaks)," peaks found in peak file(s)."))

GRanges_peaks=df_peaks %>% GRanges()

# merge peaks ----
union_peaks=c()

for(file_peak_one in files_peak){
  # print(file_peak_one)
  if(length(union_peaks) == 0){
    union_peaks=GRanges_peaks[GRanges_peaks$file==file_peak_one,]
  } else {
    union_peaks=GenomicRanges::union(union_peaks, GRanges_peaks[GRanges_peaks$file==file_peak_one,])
  }
}

message(paste0("\nTotally ",length(union_peaks)," peaks after union."))

# get counts parameters ----
blacklist = read.table(file_blacklist, sep="\t")
names(blacklist)[1:3]=c("chrom","start","end")
blacklist=blacklist %>% mutate(start = start+1)
blacklist = GRanges(blacklist)

if(species=="mouse"){standard.chr <- paste0("chr", c(1:19, "X", "Y"))} else if (species=="human"){standard.chr <- paste0("chr", c(1:22, "X", "Y"))} # only use standard chromosomes

param = readParam(max.frag=1000, pe="both", discard=blacklist, restrict=standard.chr)

# get counts parameters ----
files_bam=readLines(file_bam)
message("using csaw::regionCounts to get peak counts ...")
peak_counts = regionCounts(files_bam, union_peaks, param=param)

saveRDS(peak_counts,paste0(dir_out,"peakCounts.",label,".all.rds"))

df_row_all=rowRanges(peak_counts) %>% as.data.frame() %>% mutate(feature=paste0(seqnames,":",start,"-",end),start0=start-1)
write.table(df_row_all[,c("seqnames","start0","end")],paste0(dir_out,"peakMerged.",label,".all",".bed"),sep="\t",row.names = F,col.names = F,quote = F)

#get filtered peaks ----
cutoff=3
message(paste0("filter peaks with cutoff = ",cutoff,"\n"))
peak.abundances = edgeR::aveLogCPM(c(peak_counts))

peak_counts_filtered = peak_counts[peak.abundances > cutoff, ] # only use peaks logCPM > 3

pdf(paste0(dir_out,"histgram.peakCPM.",label,".pdf"),width=7,height = 5)
hist(peak.abundances, breaks = 5000,main="")
title(paste0(dim(peak_counts_filtered)[1]," peaks left after filtering (keep peaks with CPM >",cutoff,")"))
dev.off()

#get filtered count matrix ----
df_counts=assays(peak_counts_filtered)[["counts"]] %>% as.data.frame()

message(paste0(dim(df_counts)[1]," peaks left after filtering"))
df_info=read_tsv(file_info)

df_row=rowRanges(peak_counts_filtered) %>% as.data.frame() %>% mutate(feature=paste0(seqnames,":",start,"-",end),start0=start-1)
df_col=peak_counts_filtered@colData %>% as.data.frame() %>% mutate(file=basename(bam.files),COH_sample=word(file,1,sep="[.]")) %>% left_join(df_info)

row.names(df_counts)=df_row$feature
colnames(df_counts)=df_col$COH_sample

saveRDS(df_counts,paste0(dir_out,"peakCounts.",label,".filtered",".rds"))

write.table(df_row[,c("seqnames","start0","end")],paste0(dir_out,"peakMerged.",label,".filtered",".bed"),sep="\t",row.names = F,col.names = F,quote = F)

#get obj with count matrix ----
row.names(df_col)=df_col$COH_sample

obj_=list(SE=SummarizedExperiment(assays=list(counts=df_counts), colData = df_col))

saveRDS(obj_,paste0(dir_out,"obj.",label,".filtered",".rds"))

message("\nDone peak counts\n")











