library(Rtsne)
library(caret)
#source("C://Users/zuhu/OneDrive - City of Hope National Medical Center/statistics_R/1.functions/19.tsne/HTSeqInit.R")
file_gene_infor="C://Users/zuhu/OneDrive - City of Hope National Medical Center/statistics_R/1.functions/19.tsne/Homo_sapiens.GRCh38.V102.summary.txt"
#write_tsv ------------------------
write_tsv=function(indata,outfile){
  write.table(indata,outfile,col.names = T,row.names = F,quote = F,na = "",sep="\t")
}

#read tsv --------------
read_tsv=function(infile){
  read.table(infile,header = T,sep = "\t",stringsAsFactors = F)
}


# Creat Deseq object from dataset -----------------------------------------------------------------------
#data_diag is can be a filtered sample list. 
#this function will select the collapsed samples. Default sample name should be "COH_sample"
#var_effect can be vetor of diagnosis or/and batch

create_deseqObject_from_matrixdata=function(data_diag,data_matrix,var_effect){
  
  
  sample_in=colnames(data_matrix)
  data_diag1=data_diag[data_diag$COH_sample %in% sample_in,] %>% arrange(COH_sample)
  
  cat(paste0("Used N=",length(data_diag1$COH_sample)))
  cat("\n")
  
  data_matrix1=data_matrix[,c("feature",data_diag1$COH_sample)]
  
  #create deseq2 object
  matrix_htseq=as.matrix(data_matrix1[-1])
  row.names(matrix_htseq)=data_matrix$feature
  table(data_diag1$COH_sample==colnames(matrix_htseq))
  
  formula_in=formula(paste0("~",paste0(var_effect,collapse = "+")))
  
  htseq_object = DESeqDataSetFromMatrix(countData = matrix_htseq,
                                        colData = data_diag1,
                                        design = formula_in)
  htseq_object
  
}


# remove_high_corr_genes -----------------------------------------------------------------------


remove_high_corr_genes=function(indata,highCorrCutoff){
  set.seed(0)
  modelDF=t(indata) %>% as.data.frame()
  corrMatrix = cor(modelDF)
  highCorrIdx = findCorrelation(corrMatrix, cutoff=highCorrCutoff)
  highCorrGenes=colnames(modelDF[, highCorrIdx])
  
  modelDFnoCorr=modelDF[, -highCorrIdx]
  noHighCorrGeneDF=t(modelDFnoCorr)
  noHighCorrGeneDF
}


# draw density plot
rlogDist=function(rlogDF, outDir, fout, threshLine=0.3){
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


# generate_data_for_tsne_from_DeseqObejct_NoBatchCorrection------------------------------------
generate_data_for_tsne_from_DeseqObejct_NoBatchCorrection=function(in_object,file_gene_infor,dir_density,file_name_density){
  
  #VarianceStabilizingTransformation
  cat("Running vst\n")
  htseq_vst=DESeq2::varianceStabilizingTransformation(in_object,blind = T)
  htseq_rlog=round(data.frame(assay(htseq_vst)), digits = 2 )
  
  cat("Draw density plot\n")
  
  if(!dir.exists(dir_density)){dir.create(dir_density)}
  rlogDist(htseq_rlog, dir_density,fout = file_name_density, threshLine = 0.3)
  
  #extract coding gene
  cat("extract coding gene\n")
  info_gtf=read.table(file_gene_infor,header = T,stringsAsFactors = F)
  codinggeneList=info_gtf$gene_id[info_gtf$gene_biotype=="protein_coding"]
  htseq_rlog=htseq_rlog[row.names(htseq_rlog) %in% codinggeneList,]
  
  
  #batch effect correction
  
  #modcombat = model.matrix(~1, data=data_diag)
  #htseq_correctBatch = ComBat(dat=as.matrix(htseq_rlog), batch=data_diag$batch, mod=modcombat, par.prior=TRUE) %>%
  #  as.data.frame() %>% round(., digits = 2)
  
  
  #remove genes on chr X, Y and MT
  cat("remove genes on chr X, Y and MT\n")
  geneXYM=info_gtf$gene_id[info_gtf$chr %in% c('chrM','chrX','chrY')]
  htseq_noXYM=htseq_rlog[!row.names(htseq_rlog) %in% geneXYM,]
  
  #remove consistantly low-expressed genes
  htseq_highExpr=htseq_noXYM[apply(htseq_noXYM, 1, max) >= 5,]
  
  cat("running MAD\n")
  #remove low mad genes
  N=10000
  geneMadOrd = names(sort(apply(htseq_highExpr,1,mad),decreasing = T))
  topMadGenes=geneMadOrd[1:N]
  htseq_topMAD=htseq_highExpr[topMadGenes,]
  
  #remove high corr genes, find attributes that are highly corrected (ideally >0.75)
  cat("Removing correlated genes\n")
  indata=htseq_topMAD
  highCorrCutoff=0.75
  htseq_forPlot=remove_high_corr_genes(htseq_topMAD,0.75)
  dim(htseq_forPlot)
  htseq_forPlot
  
}




#draw tsne plot
  
run_tsne=function(indata,N,i){
  gene_names_decSort=names(sort(apply(indata,1,mad),decreasing = T))
  expMadGenes = gene_names_decSort[1:N]
  topMadDataTv=t(indata[expMadGenes,])
  cat("Perplexity:", i, "; top MAD gene number:", dim(topMadDataTv)[2], "\n");
  dim(topMadDataTv)
  # Run TSNE
  set.seed(10) # Sets seed for reproducibility
  rlogTsneOut = Rtsne(topMadDataTv, dims = 2, perplexity = i, 
                      theta = 0.1, 
                      max_iter = 5000, 
                      check_duplicates = F,
                      partial_pca=T, 
                      num_threads = 4)
  
  tsneDF=data.frame(COH_sample=row.names(topMadDataTv)) %>% 
    mutate(X=rlogTsneOut$Y[,1], 
           Y=rlogTsneOut$Y[,2],
           perplexity=i,
           gene_n=N) 
  
  tsneDF
}



draw_tsne=function(indata,data_diag,diag_col,N,i,axis_by){
  data_plot=indata %>% left_join(data_diag) %>% na.omit()
  
  ggplot() + xlab(paste0("top MAD gene num = ", N, "; perplexity=", i)) + ylab("")+
    theme_bw() +
    geom_point(data=data_plot, aes(X, Y, color=diag)) +
    scale_color_manual(values = diag_col) +
    scale_x_continuous(breaks = round(seq(floor(min(data_plot$X)), ceiling(max(data_plot$X)), by = axis_by),1)) +
    scale_y_continuous(breaks = round(seq(floor(min(data_plot$Y)), ceiling(max(data_plot$Y)), by = axis_by),1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) )
}









