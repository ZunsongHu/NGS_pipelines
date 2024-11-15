library(Rtsne)
library(caret)
library(DESeq2)
library(tidyr)
library(dplyr)
#source("C://Users/zuhu/OneDrive - City of Hope National Medical Center/statistics_R/1.functions/19.tsne/HTSeqInit.R")

#---

n_pool=c(100,150,300,600,1000,1200)
i_pool=c(10,20,30,40)

print_index=function(){
  for(N in n_list){
    for (i in i_pool){
      run_index=(match(N,n_list)-1)*length(i_pool)+match(i,i_pool)
      print(run_index)
    }}
}

generate_html_stastics=function(indata,file,var_list){
  summarytools::view(summarytools::dfSummary(indata[var_list]),file,method = "panel")
}


#write_tsv ------------------------
write_tsv=function(indata,outfile){
  write.table(indata,outfile,col.names = T,row.names = F,quote = F,na = "",sep="\t")
}

#read tsv --------------
read_tsv=function(infile){
  read.table(infile,header = T,sep = "\t",stringsAsFactors = F)
}
# gene list in --------------------------------

file_gene_infor="C://Users/zuhu/OneDrive - City of Hope National Medical Center/statistics_R/1.functions/19.tsne/Homo_sapiens.GRCh38.V102.summary.txt"
dux4_genelist=readLines("C://Users/zuhu/OneDrive - City of Hope National Medical Center/project/ZGu/B_ALL/RNAseq/3.code/1.prepare_matrix/2.data_temp/dux4_gene.list")

info_gtf=read.table(file_gene_infor,header = T,stringsAsFactors = F)
info_gtf_dux4=info_gtf[info_gtf$gene_name %in% dux4_genelist,]
codinggeneList=info_gtf$gene_id[info_gtf$gene_biotype=="protein_coding" | info_gtf$gene_name %in% dux4_genelist]
geneXYM=info_gtf$gene_id[info_gtf$chr %in% c('chrM','chrX','chrY')]
gene_in=codinggeneList[!codinggeneList %in% geneXYM]

rm(file_gene_infor)
rm(info_gtf_dux4)
rm(codinggeneList)
rm(geneXYM)



########################################  Calculate freq, no group variable, one or multiple variables ########################################

#value_x=unlist(indata["hp"])
#var_x="hp"

cal_freq_NoGroup=function(indata,var_list,digit_n){
  bind_rows(lapply(var_list, function(var_x){
    value_x=unlist(indata[var_x])
    out1=data.frame(table(value_x))
    out2=data.frame(
      Variable=var_x,
      variable_levels=paste0(as.character(unlist(out1[,1]))),
      N_used=as.character(sum(!is.na(value_x))),
      All=
        paste0(unlist(out1[,2]),
               " (",
               sprintf("%.1f",unlist(out1[,2])/sum(unlist(out1[,2]))*100),
               "%)"
        ),
      stringsAsFactors = F
    )
    out2$N_used[2:nrow(out2)]=NA
    bind_rows(data.frame(Variable=var_x,variable_levels=var_x,stringsAsFactors = F),out2) 
  }))
  
}


# Creat Deseq object from dataset -----------------------------------------------------------------------
#data_diag is can be a filtered sample list. 
#this function will select the collapsed samples. Default sample name should be "COH_sample"
#var_effect can be vetor of diagnosis or/and batch
#matrix_in must has a rowname correcsponding to the feature names!!!!!!!

create_deseqObject_from_matrixdata=function(matrix_in,diag_in,sample_var_name,design_var){
  
  sample_name_in=names(matrix_in)[names(matrix_in) %in% unlist(diag_in[sample_var_name])]
  sample_name_in=sort(sample_name_in)
  feature_in=rownames(matrix_in)
  
  cat(paste0("Used feature N=",length(feature_in)));  cat("\n")
  cat(paste0("Used sample N=",length(sample_name_in)));  cat("\n")
  
  matrix_in1=matrix_in[,sample_name_in]
  
  row.names(matrix_in1)=feature_in
  
  diag_in1=diag_in[unlist(diag_in[sample_var_name]) %in% sample_name_in,] 
  diag_in1=diag_in1[order(unlist(diag_in1[sample_var_name])),]
  
  formula_in=formula(paste0("~",paste0(design_var,collapse = "+")))
  
  htseq_object = DESeqDataSetFromMatrix(countData = matrix_in1,
                                        colData = diag_in1,
                                        design = formula_in)
  cat(dim(htseq_object));cat("\n")
  htseq_object
}


# get count matrix -------------------------
#files=list.files("1.data_raw/Counts_Ibrahim/",full.names = T)

#file=files[1]

get_count_matrix=function(files){
  for(file in files){
    id=word(basename(file),1,sep="[.]")
    count_one=read.table(file,header = F,stringsAsFactors = F) 
    names(count_one)=c("feature",id)
    if(match(file,files)==1){count_all=count_one}
    if(match(file,files)>1){count_all=count_all %>% left_join(count_one)}
  }
  count_all
}


run_vst=function(object_in){
  cat("Running vst\n")
  set.seed(0)
  object_out=DESeq2::varianceStabilizingTransformation(object_in,blind = T)
  data_round=data.frame(assay(object_out))
  dim(data_round)
  data_round
}

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
  plot(1,1, type='n', ylim=c(0,maxY), xlim=c(5,18), xlab = "Normalization value", ylab="Density", main = "" ,cex.lab=1.5,cex.axis=1.5)
  for( i in 1:sampleNum){
    points(dens.x[,i], dens.y[,i], type = "l", col=sampleCol[i]);
  }
  densX=rowMedians(dens.x)
  densY=rowMedians(dens.y)
  points(densX, densY, type='l', lwd = 2)
  #dist
  distDens=c()
  for( r in 1:sampleNum){
    distI=dist(rbind(densY, dens.y[,r]))
    distDens=c(distDens, as.vector(distI))
  }
  names(distDens)=samples
  write.table(x = distDens, file = file.path(outDir,"dist2Median.txt"), quote=F, row.names = T, sep = "\t")
  distDensSrt=sort(distDens, decreasing = F)
  plot(x = 1:sampleNum, y=distDensSrt, xlab = "Sample index", ylab="Distance to the median line",cex.lab=1.5,cex.axis=1.5)
  abline(h = threshLine)
  text(x = 5, y=threshLine, labels = paste("dist =",threshLine), pos = 3)
  dev.off()
  par(mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1)
}

draw_density=function(matrix_in,threshLine,dir_out){
  #matrix_in: col, samples. raw, features. row name of matrix_in must be feature names.
  if(!dir.exists(dir_out)){dir.create(dir_out)}
  
  df_sample=data.frame(sample=colnames(matrix_in),stringsAsFactors = F)
  df_sample$col=rainbow(nrow(df_sample), alpha = 0.6)
  
  sample_list=unique(df_sample$sample)
  
  df_density=bind_cols(lapply(sample_list, function(sample_one){
    density_one=density(unlist(matrix_in[sample_one]))
    df_density_one=data.frame(x=round(density_one$x,4),y=round(density_one$y,4))
    names(df_density_one)=paste0(sample_one,"_",names(df_density_one))
    df_density_one
  }))
  
  max_x=max(unlist(df_density[which(grepl("x",names(df_density)))]))
  min_x= min(unlist(df_density[which(grepl("x",names(df_density)))]))
  
  max_y=max(unlist(df_density[which(grepl("y",names(df_density)))]))
  min_y= min(unlist(df_density[which(grepl("y",names(df_density)))]))
  
  file_name=file.path(dir_out,paste0("density_sampleN",length(sample_list),"_featureN",nrow(matrix_in),".tsv"))
  write_tsv(df_density,file_name)
  
  pdf_name_dens=paste0("density_sampleN",length(sample_list),"_featureN",nrow(matrix_in),".pdf")
  
  pdf(file = file.path(dir_out,pdf_name_dens), width = 10, height = 6)
  par(mar=c(5,5,2,2))
  plot(NA,xlim=c(min_x,max_x),ylim=c(min_y,max_y),xlab = "Batch effect corrected normalized value", ylab="Density", main = "" ,cex.lab=1.5,cex.axis=1.5)
  
  for(sample_one in sample_list){
    x_value=unlist(df_density[paste0(sample_one,"_x")])
    y_value=unlist(df_density[paste0(sample_one,"_y")])
    line_col=df_sample$col[df_sample$sample==sample_one]
    lines(x=x_value,y=y_value,col=line_col)
  }
  dev.off()
  
  #dist plot
  df_y=df_density[names(df_density)[grepl("y",names(df_density))]]
  
  median_value=rowMedians(as.matrix(df_y))
  
  df_dist=bind_rows(lapply(sample_list, function(sample_one){
    data.frame(
      sample=sample_one,
      dist_value=as.vector(dist(t(data.frame(median_value=median_value,
                                             y_value=unlist(df_density[paste0(sample_one,"_y")]))))),
      stringsAsFactors = F
    )
  }))
  
  
  
  df_dist=df_dist %>% arrange(dist_value) %>% mutate(obs=1:nrow(df_dist)) %>% arrange(desc(dist_value))
  
  file_name_dist=file.path(dir_out,paste0("distance_sampleN",length(sample_list),"_featureN",nrow(matrix_in),".tsv"))
  #write_tsv(df_dist,file_name_dist)
  
  
  pdf_name_dist=paste0("distance_sampleN",length(sample_list),"_featureN",nrow(matrix_in),".pdf")
  #pdf(file = file.path(dir_out,pdf_name_dist), width = 10, height = 6)
  #par(mar=c(5,5,2,2))
  
  #plot(df_dist$obs,df_dist$dist_value,xlab = "Sample index", ylab="Distance to the median line", main = "" ,cex.lab=1.5,cex.axis=1.5)
  #abline(h = threshLine)
  #text(x = 5, y=threshLine, labels = paste("dist =",threshLine), pos = 3)
  #dev.off()
}


draw_density_plot=function(indata,dir_out,file_density,dist_cutoff){
  if(!dir.exists(dir_out)){dir.create(dir_out)}
  rlogDist(indata, dir_out,fout = file_density, threshLine = dist_cutoff)
  data_dist=read.table(paste0(dir_out,"/dist2Median.txt"),skip = 1,stringsAsFactors = F)
  names(data_dist)=c("sample_id","dist")
  data_dist
}


get_high_MAD_gene=function(indata,N_genes){
  geneMadOrd = names(sort(apply(indata,1,mad),decreasing = T))
  if(length(geneMadOrd)>N_genes){topMadGenes=geneMadOrd[1:N_genes]}
  if(length(geneMadOrd)<=N_genes){topMadGenes=geneMadOrd}
  outdata=indata[topMadGenes,]

}

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


remove_nonCOdingXYM_genes=function(indata){
  info_gtf=read.table(file_gene_infor,header = T,stringsAsFactors = F)
  codinggeneList=info_gtf$gene_id[info_gtf$gene_biotype=="protein_coding"]
  geneXYM=info_gtf$gene_id[info_gtf$chr %in% c('chrM','chrX','chrY')]
  gene_in=codinggeneList[!codinggeneList %in% geneXYM]
  outdata=indata[row.names(indata) %in% gene_in,]
}

remove_lowexpression_genes=function(indata,cutoff){
  out_data=indata[apply(indata, 1, max) >= cutoff,]
}

#batch effect correction -------------------------
correct_batch=function(indata,diag_in,sample_id){
  print(dim(indata))
  print("running batch correction")
  
  sample_matrix=names(indata)
  sample_diag=unlist(diag_in[sample_id])
  
  sample_in=sort(sample_matrix[sample_matrix %in% sample_diag])
  
  data_diag1=diag_in[unlist(diag_in[sample_id]) %in% sample_in,] %>% arrange(sample_id) %>%
    mutate(batch1=word(batch,1,sep="_"),batch2=word(batch,2,sep="_"),batch3=word(batch,3,sep="_"))
  
  indata1=indata[,unlist(data_diag1[sample_id])]
  
  modcombat = model.matrix(~1, data=data_diag1)
  
  print("1. Correct for mRNA/totalRNA")
  out_data1 = sva::ComBat(dat=as.matrix(indata1), batch=data_diag1$batch2, mod=modcombat, par.prior=TRUE) %>% as.data.frame() 
  print("2. Correct for stranded/non-stranded")
  out_data2 = sva::ComBat(dat=as.matrix(out_data1), batch=data_diag1$batch1, mod=modcombat, par.prior=TRUE) %>% as.data.frame() 
  print("3. Correct for length:100bp/150bp/75bp etc.")
  out_data3 = sva::ComBat(dat=as.matrix(out_data2), batch=data_diag1$batch2, mod=modcombat, par.prior=TRUE) %>% as.data.frame() %>% round(., digits = 5)
  
  print(paste0("outdata dim: ", paste0(dim(out_data3),collapse = ",")))
  out_data3
}




# for tsne ----------------------
from_vst2tsnePanel_df=function(df_in,diag_in){
  
  df_vst_batchC=correct_batch(df_in,diag_in,"COH_sample","batch")
  
  df_vst_batchC_codingNoXYM=df_vst_batchC[row.names(df_vst_batchC) %in% gene_in,]
  
  df_vst_batchC_codingNoXYM_highE=remove_lowexpression_genes(df_vst_batchC_codingNoXYM,5)
  
  df_vst_batchC_codingNoXYM_highE_MAD3000=get_high_MAD_gene(df_vst_batchC_codingNoXYM_highE,3000)
  
  df_vst_batchC_codingNoXYM_highE_MAD3000_noCorr=remove_high_corr_genes(df_vst_batchC_codingNoXYM_highE_MAD3000,0.75)
  
  df_vst_batchC_codingNoXYM_highE_MAD3000_noCorr_MAD1000=get_high_MAD_gene(df_vst_batchC_codingNoXYM_highE_MAD3000_noCorr,1500)
  print(dim(df_vst_batchC_codingNoXYM_highE_MAD3000_noCorr_MAD1000))
  df_vst_batchC_codingNoXYM_highE_MAD3000_noCorr_MAD1000
}

from_vstbatch2tsnePanel_df=function(df_in,gene_in){
  print("codingNoXYM")
  df_vst_batchC_codingNoXYM=df_in[row.names(df_in) %in% gene_in,]
  print("remove low exp")
  df_vst_batchC_codingNoXYM_highE=remove_lowexpression_genes(df_vst_batchC_codingNoXYM,5)
  print("top MAD10000")
  df_vst_batchC_codingNoXYM_highE_MAD10000=get_high_MAD_gene(df_vst_batchC_codingNoXYM_highE,10000)
  print("remove corr")
  df_vst_batchC_codingNoXYM_highE_MAD10000_noCorr=remove_high_corr_genes(df_vst_batchC_codingNoXYM_highE_MAD10000,0.75)
  print("top MAD1500")
  df_vst_batchC_codingNoXYM_highE_MAD10000_noCorr_MAD1000=get_high_MAD_gene(df_vst_batchC_codingNoXYM_highE_MAD10000_noCorr,1500)
  print(dim(df_vst_batchC_codingNoXYM_highE_MAD10000_noCorr_MAD1000))
  df_vst_batchC_codingNoXYM_highE_MAD10000_noCorr_MAD1000
}

#draw tsne plot
  
run_tsne=function(indata,i){
  indata1=t(indata)
  cat("Perplexity:", i, "; top gene number:", dim(indata1)[2], "\n");
  # Run TSNE
  set.seed(10) # Sets seed for reproducibility
  rlogTsneOut = Rtsne(indata1, dims = 2, perplexity = i, 
                      theta = 0.1, 
                      max_iter = 5000, 
                      check_duplicates = F,
                      partial_pca=T, 
                      num_threads = 4)
  
  tsneDF=data.frame(COH_sample=row.names(indata1)) %>% 
    mutate(X=rlogTsneOut$Y[,1], 
           Y=rlogTsneOut$Y[,2],
           perplexity=i,
           gene_n=dim(indata1)[2]) 
  
  tsneDF
}

run_tsne_MultipleD=function(indata,i,D){
  indata1=t(indata)
  cat("Perplexity:", i, "; top gene number:", dim(indata1)[2], "\n");
  # Run TSNE
  set.seed(10) # Sets seed for reproducibility
  rlogTsneOut = Rtsne(indata1, dims = D, perplexity = i, 
                      theta = 0.1, 
                      max_iter = 5000, 
                      check_duplicates = F,
                      partial_pca=T, 
                      num_threads = 4)
  
  tsneDF=data.frame(COH_sample=row.names(indata1)) %>% 
    mutate(X=rlogTsneOut$Y[,1], 
           Y=rlogTsneOut$Y[,2],
           Z=rlogTsneOut$Y[,3],
           perplexity=i,
           gene_n=dim(indata1)[2]) 
  
  tsneDF
}

run_tsne_panel_basedOnTopGeneList=function(indata,genelist,n_pool,i_pool){
  df_tsne_panel=bind_rows(lapply(n_pool, function(n){
    lapply(i_pool, function(i){
      feature_in=head(genelist,n)
      indata1=indata[feature_in,]
      df_tsne_temp=run_tsne_MultipleD(indata1,i,3)
      df_tsne_temp
    })
  }))
  df_tsne_panel
}


run_tsne_panel_MAD=function(indata,n_pool,i_pool){
  df_tsne_panel=bind_rows(lapply(n_pool, function(n){
    lapply(i_pool, function(i){
      indata1=get_high_MAD_gene(indata,n)
      df_tsne_temp=run_tsne(indata1,i)
      df_tsne_temp
    })
  }))
  df_tsne_panel
}


top_MAD=function(indata,n_pool){
  df_=bind_rows(lapply(n_pool, function(n){
    indata1=get_high_MAD_gene(indata,n)
    out_one=data.frame(
      feature=row.names(indata1),
      group=paste0("Top",n),
      stringsAsFactors = F
    )
    out_one
  }))
  df_
}



run_tsne_panel_boruta=function(indata,boruta_imp,n_pool,i_pool){
  boruta_imp1=boruta_imp %>% arrange(desc(meanImp)) %>% filter(decision=="Confirmed")
  df_tsne_panel=bind_rows(lapply(n_pool, function(n){
    lapply(i_pool, function(i){
      boruta_imp2=boruta_imp1[1:n,]
      indata1=indata[row.names(indata) %in% row.names(boruta_imp2),]
      df_tsne_temp=run_tsne(indata1,i)
      df_tsne_temp
    })
  }))
  df_tsne_panel
}

run_tsne_boruta=function(indata,N,i){
  gene_in=row.names(importance_boruta_in)[1:N]
  
  data_for_tsne=t(indata[gene_in,])
  
  
  cat("Perplexity:", i, "; top MAD gene number:", dim(data_for_tsne)[2], "\n");
  dim(data_for_tsne)
  # Run TSNE
  set.seed(10) # Sets seed for reproducibility
  rlogTsneOut = Rtsne(data_for_tsne, dims = 2, perplexity = i, 
                      theta = 0.1, 
                      max_iter = 5000, 
                      check_duplicates = F,
                      partial_pca=T, 
                      num_threads = 4)
  
  tsneDF=data.frame(COH_sample=row.names(data_for_tsne)) %>% 
    mutate(X=rlogTsneOut$Y[,1], 
           Y=rlogTsneOut$Y[,2],
           perplexity=i,
           gene_n=N) 
  
  tsneDF
}

draw_tsne_MAD=function(indata,data_diag,diag_col,N,i,axis_by){
  data_plot=indata %>% left_join(data_diag) %>% drop_na()
  
  ggplot() + 
    xlab(paste0("tSNE dimension 1\ntop gene num = ", N, "; perplexity=", i)) +
    ylab("tSNE dimension 2")+
    theme_bw() +
    geom_point(data=data_plot, aes(X, Y, color=diag)) +
    scale_color_manual(values = diag_col) +
    scale_x_continuous(breaks = round(seq(floor(min(data_plot$X)), ceiling(max(data_plot$X)), by = axis_by),1)) +
    scale_y_continuous(breaks = round(seq(floor(min(data_plot$Y)), ceiling(max(data_plot$Y)), by = axis_by),1)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15,face="bold")
    ) +
    scale_fill_discrete(name="Experimental\nCondition")
}


draw_tsne_MAD_panel=function(indata,i_pool,file_path,file_name,col.n=6,legend_loc="bottom"){
  n_list=unique(indata$gene_n)
  n_plot=length(n_list)*length(i_pool)
  print(paste0("plot number:",n_plot))
  col.r=n_plot/col.n
  
  tsnePanelList=list()
  for (i in i_pool){
    for(N in n_list){
      run_index=(match(N,n_list)-1)*length(i_pool)+match(i,i_pool)
      print(run_index)
      df_tsne_in=indata %>% filter(gene_n==N & perplexity==i)
      tsnePanelList[[run_index]]=draw_tsne_MAD(df_tsne_in,diag_in_tsne,subtypeCol,N,i,10)
    }}
  
  tsnePanel=ggarrange(plotlist=tsnePanelList,ncol=col.n,nrow=col.r,common.legend = T,legend=legend_loc)
  file_=file.path(file_path,paste0("tSNE_panel_",file_name))
  ggsave(paste0(file_,".pdf"),width = 5*col.n,height = 4.5*col.r)
  ggsave(paste0(file_,".png"),width = 5*col.n,height = 4.5*col.r,dpi = 300)
}



# color -------------------
#BALL color
subtypeCol=c()
{
  subtypeCol["ETV6-RUNX1"]="gold2"
  subtypeCol["ETV6-RUNX1-like"]="pink"
  subtypeCol["KMT2A"]="#1F78B5"
  subtypeCol["Ph"]="magenta3"
  subtypeCol["DUX4"]='grey40'
  subtypeCol["TCF3-PBX1"]="darkgoldenrod4"
  subtypeCol["ZNF384"]="#A8DD00"
  subtypeCol["MEF2D"]="#66C2A6"
  subtypeCol["BCL2/MYC"]="seagreen2"
  subtypeCol["NUTM1"]='black'
  subtypeCol["HLF"]= "skyblue"
  subtypeCol["PAX5(P80R)"]="orangered"
  subtypeCol["PAX5 P80R"]="orangered"
  subtypeCol["Hyperdiploid"]="#3E9F32"
  subtypeCol["LowHypo"]="#1E90FF"
  subtypeCol["Low hypodiploid"]="#1E90FF"
  subtypeCol["NearHaploid"]='blue3'
  subtypeCol["Near haploid"]='blue3'
  subtypeCol["Ph-like"]="red4"
  subtypeCol["PAX5alt"]="#FFA620"
  subtypeCol["PAX5-ETV6"]="#808000"
  subtypeCol["iAMP21"]="lightslateblue"
  subtypeCol["IKZF1(N159Y)"]="#CCCC33"
  subtypeCol["IKZF1 N159Y"]="#CCCC33"
  subtypeCol["LowHyper"]="cyan"
  subtypeCol["Bother"]='grey75'
  subtypeCol["Low hyperdiploid"]='grey75'
  subtypeCol["CRLF2(non-Ph-like)"]='grey75'
  subtypeCol["KMT2A-like"]='grey75'
  subtypeCol["ZNF384-like"]='grey75'
  subtypeCol["Other"]='grey75'
  subtypeCol["ZEB2/CEBPE"]="#D27B1C86"
  subtypeCol["Y"]="#E6BEFF"
  subtypeCol["Unknown"]="#469990"
  subtypeCol["_Prediction"]="red4"
}


df_ball_col=data.frame(
  diag=names(subtypeCol),
  col=subtypeCol,
  obs=1:length(subtypeCol),
  stringsAsFactors = F
)




