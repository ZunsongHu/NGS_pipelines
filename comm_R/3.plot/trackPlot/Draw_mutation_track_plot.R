suppressMessages(library(Gviz))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))

#temp parameters ------------------

fileBam="/scratch/zuhu/project/zhaohuigu/singlecell/R_code/mutation/COH002922_newCounnt_total.CB.cutoff10/mutation_bam/NRAS.G12S/AACGTTGGTCTGATTG-1.bam"
dir_out="/scratch/zuhu/project/zhaohuigu/singlecell/R_code/mutation/COH002922_newCounnt_total.CB.cutoff10/mutation_bam/NRAS.G12S/"

#args <- commandArgs(trailingOnly = TRUE)
#print(args)

#file_tsv=args[1]

options(ucscChromosomeNames=FALSE)


gene="NRAS"
mut="G13D"
mutation_position=114716127

interval=50
chr="chr1"
region_start=mutation_position-interval
region_end=mutation_position+interval


genomeAxis = GenomeAxisTrack(name=paste0(gene,".",mut)) 
sTrack = SequenceTrack(BSgenome.Hsapiens.UCSC.hg38) 
customFromTxDb = GeneRegionTrack(TxDb.Hsapiens.UCSC.hg38.knownGene,chromosome=chr,from=region_start,to=region_end) 


activatedReads = AlignmentsTrack(fileBam, isPaired = TRUE, referenceSequence=sTrack) 

png(paste0(dir_out,"tmp.png"),height=2500,width=2500,res=300)
plotTracks(c(activatedReads,sTrack,genomeAxis),from=region_start,to=region_end,labelPos="below",min.height = 8,cex=0.5) 
dev.off()


pdf(paste0(dir_out,"tmp.pdf"),height=10,width=10)
plotTracks(c(activatedReads,sTrack,genomeAxis),from=region_start,to=region_end,labelPos="below",min.height = 8,cex=0.5) 
dev.off()




















