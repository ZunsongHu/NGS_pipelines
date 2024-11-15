#read in reference data
library(DESeq2)
library(tidyr)
library(dplyr)
library(readr)
library(magrittr)
#refDir='./etc/';
#GTF=read_tsv(file.path(refDir,'Homo_sapiens.GRCh37.V75.summary.txt') )
#id2gene=GTF %>% select(ID=EGID, gene=Symbol)
#codingGenes=GTF %>% filter(Type == "protein_coding") %>% select(ID=EGID, gene=Symbol)
#chrXYM_Genes=GTF %>% filter(Chr == "X" | Chr == "Y" | Chr == "MT")

#convert HTSeq to count and rld
#keep the genes with at least N samples covered by >= 10 reads
HTSeqConvert <- function(sampleInfor, N=1){
  print("Convert HTSeq to pass QC genes!");
  sampleTable=data.frame(sampleName=sampleInfor$sample,
                         fileName=sampleInfor$path)
  ddsHTSeq_0=DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~ 1)
  print("Loading HTSeq files is done!")
  
  #keep the genes with at least N samples covered by >= 10 reads
  ddsHTSeq <- ddsHTSeq_0[ apply(X = counts(ddsHTSeq_0)
                                , MARGIN = 1
                                , FUN = function(x) sum(x>=10) ) >=N, ]
  rldHTSeq <- varianceStabilizingTransformation(ddsHTSeq, blind=T)
  print("VST transformation is done!")
  return(list(count=counts(ddsHTSeq), rld=rldHTSeq, dds=ddsHTSeq))
}

#convert HTSeq to count and rld
HTSeqConvertCodingGene <- function(sampleInfor){
  print("Convert HTSeq to coding genes!");
  sampleTable=data.frame(sampleName=sampleInfor$sample,
                         fileName=sampleInfor$path)
  ddsHTSeq_0=DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~ 1)
  print("Loading HTSeq files is done!")
  
  geneNames=rownames(ddsHTSeq_0)
  ddsHTSeq=ddsHTSeq_0[which(geneNames %in% codingGenes$ID), ]
  ddsCount=counts(ddsHTSeq)
  rldHTSeq=varianceStabilizingTransformation(ddsHTSeq, blind=T)
  print("VST transformation is done!")
  return(list(count=ddsCount, rld=rldHTSeq, dds=ddsHTSeq))
}

#convert HTSeq to count and rld and keep all the genes
HTSeqConvertAllGene <- function(sampleInfor){
  print("Convert HTSeq to all genes!");
  sampleTable=data.frame(sampleName=sampleInfor$sample,
                         fileName=sampleInfor$path)
  ddsHTSeq=DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, design= ~ 1)
  ddsCount <- counts(ddsHTSeq)
  print("Loading HTSeq files is done!")
  rldHTSeq <- varianceStabilizingTransformation(ddsHTSeq, blind=T)
  print("VST transformation is done!")
  return(list(count=ddsCount, rld=rldHTSeq, dds=ddsHTSeq))
}

rlogDist <- function(rlogDF, outDir, fout, threshLine=0.3){
  samples=colnames(rlogDF)
  sampleNum=length(samples);
  sampleCol=rainbow(sampleNum, alpha = 0.6)
  dens.x=c()
  dens.y=c()
  maxY=0;
  for( i in 1:sampleNum){
    sampleI=rlogDF[,i]
    density=density(sampleI, from=-5, to=20 )
    maxY=ifelse(maxY > max(density$y), maxY, max(density$y));
    dens.x=cbind(dens.x, density$x)
    dens.y=cbind(dens.y, density$y)
  }
  pdf(file = file.path(outDir,fout), width = 10, height = 12)
  par(mfrow=c(2,1), mar=c(5,5,2,2))
  plot(1,1, type='n', ylim=c(0,1.2*maxY), xlim=c(-2,18), xlab = "rlog value", ylab="Density", main = "" )
  for( i in 1:sampleNum){
    points(dens.x[,i], dens.y[,i], type = "l", col=sampleCol[i]);
  }
  densX=rowMedians(dens.x)
  densY=rowMedians(dens.y)
  points(densX, densY, type='l', lwd = 2)
  ###########################################################
  distDens=c()
  for( r in 1:sampleNum){
    distI=dist(rbind(densY, dens.y[,r]))
    distDens=c(distDens, as.vector(distI))
  }
  names(distDens)=samples
  write.table(x = distDens, file = file.path(outDir,"dist2Median.txt"), quote=F, row.names = T, sep = "\t")
  distDensSrt=sort(distDens, decreasing = F)
  plot(x = 1:sampleNum, y=distDensSrt, xlab = "Sample index", ylab="Distance to the median line")
  abline(h = threshLine)
  text(x = 5, y=threshLine, labels = paste("dist =",threshLine), pos = 3)
  dev.off()
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}