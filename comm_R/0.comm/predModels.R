######################################## Calculate weighted score ##################################

cal_score_weighted=function(indata,var_multiple){
  out_set=bind_cols(
    lapply(var_multiple, function(var_one){
      varlue=unlist(indata[var_one])
      beta=rr_train$beta[rr_train$varname==var_one]
      value_new=varlue*beta
      out=data.frame(value_new=value_new,stringsAsFactors = F)
      names(out)=paste0(var_one,"_weighted")
      out
    }
    ))
  out_set$ziduan=indata$ziduan
  out_set$score_weighted=apply(out_set[paste0(var_multiple,"_weighted")],1,sum)
  out_set
}




#covert_betase_digitn  ---------------------------------
covert_betase_digitn=function(x,digit_n){
  
  x1=as.numeric(word(x,1,sep=" "))
  x2=as.numeric(gsub("[()]","",word(x,2,sep=" ")))
  
  x11=sprintf(paste0("%.",digit_n,"f"),x1)
  x21=sprintf(paste0("%.",digit_n,"f"),x2)
  paste0(x11," (",x21,")")
}

