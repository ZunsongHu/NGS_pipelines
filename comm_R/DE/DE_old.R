
# DE analysis, one group vs all the others ----------------------
#DE analysis is based on raw data, no normalization is needed

DE_oneGroup_vs_Others=function(in_count_matrix,in_diag,test_level){
  in_diag1=in_diag %>% mutate(diag1=ifelse(diag==test_level,diag,"0other"))
  in_diag1$diag1=as.factor(in_diag1$diag1)
  
  deseq_object=create_deseqObject_from_matrixdata(in_count_matrix,in_diag1,"COH_sample","diag1")
  dim(deseq_object)
  
  dds = DESeq2::DESeq(deseq_object)
  res=results(dds, contrast=c("diag1",test_level,"0other"))
  
  res_DE=as.data.frame(res@listData) %>% 
    mutate(feature=res@rownames,
           test=paste0(test_level,"_vs_","other")) %>% 
    arrange(padj) 
  res_DE
}

