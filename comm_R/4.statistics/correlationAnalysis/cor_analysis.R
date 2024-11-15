######################################## draw correlation matrix plot ########################################
#indata=data_for_cor
#cor_method="spearman"

draw_cor_matrix=function(indata,cor_method){
  M<-cor(indata,method = cor_method,use = "complete.obs")
  cor.mtest <- function(mat) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        tmp <- cor.test(mat[, i], mat[, j], method = cor_method)
        p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
      }
    }
    colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
    p.mat
  }
  # matrix of the p-value of the correlation
  p.mat <- cor.mtest(indata)
  
  col <- colorRampPalette(c("#4477AA","#77AADD","grey90", "#EE9988","#BB4444"  ))
  corrplot::corrplot(M, method="color", col=col(200),  
                     type="upper", order="original", 
                     addCoef.col = "black", # Add coefficient of correlation
                     tl.col="black", tl.srt=45, #Text label color and rotation
                     # Combine with significance
                     p.mat = p.mat, sig.level = 0.05, insig = "blank", 
                     # hide correlation coefficient on the principal diagonal
                     diag=FALSE 
  )

}



######################################## Remove highly correlated variables ########################################




remove_highly_corelated_var=function(indata,var_id,method,cutoff){
  
  value_data=indata[names(indata)[-match(var_id,names(indata))]]
  
  cor_out_raw=cor(value_data,method = method)
  diag(cor_out_raw)=NA
  
  if(any(abs(cor_out_raw)>=cutoff,na.rm = T)){
    out_cor_all=data.frame(
      name1=NA,
      name2=NA,
      cor_value=NA,
      stringsAsFactors = F
    )
    
    for(i in 1:nrow(cor_out_raw)){
      for( j in 1:ncol(cor_out_raw)){
        cor_one=cor_out_raw[i,j]
        if(!is.na(cor_one) & abs(cor_one)>=cutoff){
          name1=row.names(cor_out_raw)[i]
          name2=colnames(cor_out_raw)[j]
          
          name12=c(name1,name2)[order(c(name1,name2))]
          
          out_one=data.frame(
            name1=name12[1],
            name2=name12[2],
            cor_value=round(cor_one,6),
            stringsAsFactors = F
          )
          out_cor_all=rbind(out_cor_all,out_one)
        }
      }
    }
  }
  
  out_cor_all=out_cor_all %>% distinct() %>% filter(!is.na(name1))
  
  drop_var_list=unique(out_cor_all$name2)
  
  data_out=indata[names(indata)[-match(drop_var_list,names(indata))]]
  data_out
}











