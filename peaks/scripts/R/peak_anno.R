# set ---------------------------
# library(GenomicRanges)
library(ggplot2)

source("/home/zgu_labs/bin/R/functions.R")
source("/home/zuhu/bin/R/plots/F_plots.R")

# testing parameters ----
setwd("/scratch/zuhu/project/Grp_NadiaCarlesso/atac/")

dir_out=     "test/"
id=          "test"
file_peak=   "out_raw/PBS_DN_24hr/ATAC-seq/6.broadPeak/PBS_DN_24hr_peaks.filt.broadPeak"

tss="out_raw/PBS_DN_24hr/ATAC-seq/6.broadPeak/PBS_DN_24hr_peaks.filt.broadPeak.tss.bed"
ctcf="out_raw/PBS_DN_24hr/ATAC-seq/6.broadPeak/PBS_DN_24hr_peaks.filt.broadPeak.ctcf.bed"
enhancer="out_raw/PBS_DN_24hr/ATAC-seq/6.broadPeak/PBS_DN_24hr_peaks.filt.broadPeak.enhancer.bed"
promoter="out_raw/PBS_DN_24hr/ATAC-seq/6.broadPeak/PBS_DN_24hr_peaks.filt.broadPeak.promoter.bed"

homer="out_raw/PBS_DN_108hr/ATAC-seq/6.broadPeak/PBS_DN_108hr_peaks.filt.broadPeak.anno.tsv"

species_=     "mouse"

# parameters ----
args <- commandArgs(trailingOnly = TRUE)
print(args)

dir_out=      args[1]
id=           args[2]
file_peak=    args[3]

tss=          args[4]
ctcf=         args[5]
enhancer=     args[6]
promoter=     args[7]

homer=        args[8]

species_=     args[9]

cutoff=1000
# create dir out ----
if(!dir.exists(dir_out)){dir.create(dir_out)}

# read files ----
df_peak=vroom::vroom(file_peak,col_names = F)

df_tss=vroom::vroom(tss,col_names = F)
df_ctcf=vroom::vroom(ctcf,col_names = F)
df_enhancer=vroom::vroom(enhancer,col_names = F)
df_promoter=vroom::vroom(promoter,col_names = F)

#get peaks annotation ----
get_anno_peak=function(df_peak,df_one,label){
  df_one$d=unlist(df_one[paste0("X",ncol(df_one))])
  
  df_one1=df_one %>% 
    mutate(feature=paste0(X1,":",X2,"-",X3)) %>% dplyr::select(feature,d)
  
  df_peak1=df_peak %>% 
    transmute(feature=paste0(X1,":",X2,"-",X3)) %>% 
    left_join(df_one1) %>% 
    group_by(feature) %>% arrange(desc(d)) %>% 
    slice_head(n=1) %>% ungroup() %>% 
    mutate(inRegion=ifelse(d<= cutoff,'Yes',"no"))
  
  names(df_peak1)[2:3]=c(paste0("distance2",label),paste0("Near_",label))
  df_peak1
}

df_anno_tss=get_anno_peak(df_peak,df_tss,"TSS")
df_anno_ctcf=get_anno_peak(df_peak,df_ctcf,"CTCF")
df_anno_enhancer=get_anno_peak(df_peak,df_enhancer,"Enhancer")
df_anno_Promoter=get_anno_peak(df_peak,df_promoter,"Promoter")

# table(df_anno_tss$Near_TSS)
# table(df_anno_ctcf$Near_CTCF)
# table(df_anno_enhancer$Near_Enhancer)
# table(df_anno_Promoter$Near_Promoter)

df_anno_=df_anno_tss %>% left_join(df_anno_ctcf) %>% left_join(df_anno_enhancer) %>% left_join(df_anno_Promoter)

#output ----
write_tsv(df_anno_,paste0(dir_out,id,".peak.anno.tsv"))

#get anno statistics ----
df_anno_1=df_anno_ 
df_anno_1[c( "distance2TSS","distance2CTCF","distance2Enhancer", "distance2Promoter")]=NULL

df_anno_1=df_anno_1%>% reshape2::melt(id="feature") %>% 
  group_by(variable) %>% mutate(n_g=n()) %>% 
  group_by(variable,value) %>% mutate(n_g_y=n()) %>% 
  mutate(feature=NULL) %>% distinct() %>% 
  mutate(per=n_g_y/n_g) %>% 
  filter(value=="Yes")

write_tsv(df_anno_1 %>% mutate(COH_sample=id),paste0(dir_out,id,".peak.anno.percentage.tsv"))

#get barplot ----
ggplot(data=df_anno_1, aes(x=variable, y=per,fill=variable)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=round(per,3)), vjust=1.6, color="white", size=3.5) +
  labs(x="Group",y="Percentage",title="Within 1kb")+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color="grey77"),
        panel.grid.minor.y = element_line(color="grey77"),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(paste0(dir_out,"barPlot.peak.anno.percentage.",id,".png"),width=6.5,height = 4)
ggsave(paste0(dir_out,"barPlot.peak.anno.percentage.",id,".pdf"),width=6.5,height = 4)

#get homer anno ----
df_homer=vroom::vroom(homer,col_names = T) %>% 
  mutate(
    feature=paste0(Chr,":",Start-1,"-",End),
    anno_type=trimws(word(Annotation,1,sep="[()]")),
    n_g=n()
  ) %>% arrange(feature) %>% 
  group_by(anno_type) %>% 
  mutate(
    COH_sample=id,
    n_g_y=n(),
         per=round(n_g_y/n_g,4)) %>% 
  dplyr::select(COH_sample,anno_type,n_g,n_g_y,per) %>% 
  distinct()
sum(df_homer$per)
write_tsv(df_homer,paste0(dir_out,id,".peak.annoHomer.percentage.tsv"))

#get homer anno pie chart ----
col_in=unique(subtypeCol)[1:nrow(df_homer)]
names(col_in)=df_homer$anno_type

df_homer1=df_homer %>% mutate(label=paste0(anno_type,"(",sprintf("%.1f",per*100),"%)")) %>% arrange(label)

ggplot(df_homer1, aes(x="", y=per, fill=anno_type)) +
  geom_bar(stat="identity", width=1, color="white") +
  scale_fill_manual(values=col_in,labels = df_homer1$label,name="HomerAnnoGroup") +
  coord_polar("y", start=0) +
  labs(title=id)+
  theme_void() 

ggsave(paste0(dir_out,"piePlot.peak.annoHomer.percentage.",id,".png"),width=6,height = 4)
ggsave(paste0(dir_out,"piePlot.peak.annoHomer.percentage.",id,".pdf"),width=6,height = 4)


# homers=system("ls out_raw/*/ATAC-seq/6.broadPeak/*.filt.broadPeak.anno.tsv",intern = T)
# homer=homers[1]
# for (homer in homers){
#   print(match(homer,homers))
#   id=gsub("_peaks","",word(basename(homer),1,sep="[.]"))
#   dir_out=paste0(dirname(homer),"/")
#   
#   df_homer=vroom::vroom(homer,col_names = T) %>% 
#     mutate(
#       feature=paste0(Chr,":",Start-1,"-",End),
#       anno_type=trimws(word(Annotation,1,sep="[()]")),
#       n_g=n()
#     ) %>% arrange(feature) %>% 
#     group_by(anno_type) %>% 
#     mutate(
#       COH_sample=id,
#       n_g_y=n(),
#       per=round(n_g_y/n_g,4)) %>% 
#     dplyr::select(COH_sample,anno_type,n_g,n_g_y,per) %>% 
#     distinct()
#   
#   write_tsv(df_homer,paste0(dir_out,id,".peak.annoHomer.percentage.tsv"))
#   
#   col_in=unique(subtypeCol)[1:nrow(df_homer)]
#   names(col_in)=df_homer$anno_type
#   
#   df_homer1=df_homer %>% mutate(label=paste0(anno_type,"(",sprintf("%.1f",per*100),"%)")) %>% arrange(label)
#   
#   ggplot(df_homer1, aes(x="", y=per, fill=anno_type)) +
#     geom_bar(stat="identity", width=1, color="white") +
#     scale_fill_manual(values=col_in,labels = df_homer1$label,name="HomerAnnoGroup") +
#     coord_polar("y", start=0) +
#     labs(title=id)+
#     theme_void() 
#   
#   ggsave(paste0(dir_out,"piePlot.peak.annoHomer.percentage.",id,".png"),width=6,height = 4)
#   ggsave(paste0(dir_out,"piePlot.peak.annoHomer.percentage.",id,".pdf"),width=6,height = 4)
# }








