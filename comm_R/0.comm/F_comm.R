options(warn=-1)
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})
options(warn=0)

# tsv ------------------------
write_tsv=function(indata,outfile){
  if(!dir.exists(dirname(outfile))){dir.create(dirname(outfile))}
  write.table(indata,outfile,col.names = T,row.names = F,quote = F,na = "",sep="\t")
}

read_tsv=function(infile){
  read.table(infile,header = T,sep = "\t",stringsAsFactors = F)
}

# p values ---------------------------------
get_p_sig=function(x){
  ifelse(is.na(x),"ns",ifelse(x<0.001,"***",ifelse(x<0.01,"**",ifelse(x<0.05,"*",""))))
}

convert_p=function(P_raw){
  P=ifelse(P_raw<0.001,"<0.001",
           ifelse(P_raw<0.01,sprintf("%.3f",P_raw),
                  ifelse(P_raw<0.045 | P_raw==0.045,sprintf("%.2f",P_raw),
                         ifelse(P_raw<0.055 | P_raw==0.055,sprintf("%.3f",P_raw),
                                ifelse(P_raw>0.055,sprintf("%.2f",P_raw),NA
                                )))))
  P
}




































