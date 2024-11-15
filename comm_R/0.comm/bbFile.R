

#serach bb file to extract AF with position information, one record --------------------------------------

search_bb_one=function(chr,start,end,bb_in,bigBedToBed){
  start=as.numeric(start)-1
  command=paste0(bigBedToBed," -chrom=",chr," -start=",start," -end=",end," ",bb_in," stdout")
  search_out=system(command,intern=T)
  data_search_out=data.frame(search_out=search_out,stringsAsFactors = F)
  
  get_max_maf=function(x){
    if(x ==""){x1=NA}
    if(x!=""){x1=max(as.numeric(unlist(strsplit(x,","))))}
    x1
  }
  
  data_search_out=data_search_out %>% 
    dplyr::mutate(chr=stringr::word(search_out,1,sep="\t"),
                  pos=as.numeric(stringr::word(search_out,2,sep="\t"))+1,
                  ref=stringr::word(search_out,5,sep="\t"),
                  alt=stringr::word(search_out,7,sep="\t"),
                  V10=stringr::word(search_out,10,sep="\t"),
    )
  
  data_search_out=tidyr::separate_rows(data_search_out,alt,sep=",") %>% filter(!alt=="")
  data_search_out$AF=unlist(lapply(data_search_out$V10, get_max_maf)) 
  data_search_out = data_search_out %>% dplyr::mutate(start=pos,
                                                      length_ref=nchar(ref),
                                                      length_alt=nchar(alt),
                                                      end=pos-1+ifelse(length_ref>=length_alt,length_ref,length_alt),
                                                      length_ref=NULL,length_alt=NULL,V10=NULL,search_out=NULL,pos=NULL)
  
  if(nrow(data_search_out)==0){data_search_out=data.frame(chr=chr,ref="",alt="",AF=NA,start=NA,end=NA,stringsAsFactors = F)}
  if(nrow(data_search_out)>=1){data_search_out=data_search_out}
  
  names(data_search_out)=paste0(names(data_search_out),"_database")
  data_search_out$chr=chr
  data_search_out$start=start+1
  data_search_out$end=end
  data_search_out
  
  
  
}

#serach bb file to extract AF with position information, using data set, with chr, start and end --------------------------------------
search_bb_dataset=function(indata,bb_in,bigBedToBed){
  search_out=bind_rows(lapply(1:nrow(indata), function(i){
    data_one=indata[i,]
    
    chr=data_one$chr
    start=data_one$start
    end=data_one$end
    
    search_bb_one(chr,start,end,bb_in,bigBedToBed)
  }))
  
  search_out1=indata %>% left_join(search_out) 
  search_out1
}

