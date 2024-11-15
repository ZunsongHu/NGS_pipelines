#serach vcf file to extract AF with position information, one record --------------------------------------
#chr="chr1"
#start=10489025
#end=10499925
#vcf_file="/home/zuhu/snp142Common.vcf.gz"
#bcftools="/home/zuhu/miniconda3/bin/bcftools"

search_vcf_onerecord=function(chr,start,end,vcf_file,bcftools){
  para_f="'%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%FILTER\\t%INFO/AF\\n'"
  search_region=paste0(chr,":",start,"-",end)
  
  command=paste0(bcftools," query -f ",para_f," ",vcf_file," -r ",search_region," 2> /dev/null")
  search_out=system(command,intern=T)
  
  data_search_out=data.frame(raw_out=search_out,stringsAsFactors = F)
  
  data_search_out=tidyr::separate(data_search_out, raw_out, sep = "\t",into = c("chr","pos","ref","alt","qual","filter","AF"))
  
  data_search_out=tidyr::separate_rows(data_search_out,alt,sep=",") %>% filter(!alt=="")
  
  data_search_out = data_search_out %>% dplyr::mutate(start=as.numeric(pos),
                                                      length_ref=nchar(ref),
                                                      length_alt=nchar(alt),
                                                      end=start+ifelse(length_ref>=length_alt,length_ref,length_alt)-1,
                                                      length_ref=NULL,length_alt=NULL,pos=NULL)
  
  if(nrow(data_search_out)==0){data_search_out=data.frame(chr=NA,ref=NA,alt=NA,qual=NA,filter=NA,AF=NA,start=NA,end=NA,stringsAsFactors = F)}
  if(nrow(data_search_out)>=1){data_search_out=data_search_out}
  
  names(data_search_out)=paste0(names(data_search_out),"_database")
  data_search_out$chr=chr
  data_search_out$start=start
  data_search_out$end=end
  data_search_out
}

#serach vcf file to extract AF with position information, using data set, with chr, start and end --------------------------------------

search_vcf_dataset=function(indata,vcf_file,bcftools){
  search_out=bind_rows(lapply(1:nrow(indata), function(i){
    data_one=indata[i,]
    
    chr=data_one$chr
    start=data_one$start
    end=data_one$end
    
    search_vcf_onerecord(chr,start,end,vcf_file,bcftools)
  }))
  
  search_out1=indata %>% left_join(search_out) 
  search_out1
}