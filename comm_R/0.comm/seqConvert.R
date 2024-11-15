
# Convert Amino acid --------------------------------

convert_aa <- function(aa_vec) {
  # Create a lookup table of amino acid abbreviations
  aa_table <- data.frame(
    aa_3 = c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
             "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val","*"),
    aa_1 = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V","*"),
    stringsAsFactors = FALSE
  )
  
  # Use the lookup table to convert the 3-letter abbreviations to 1-letter abbreviations
  aa_vec_1 <- aa_table$aa_1[match(aa_vec, aa_table$aa_3)]
  
  return(aa_vec_1)
}

convert_aaPos=function(x){
  x=trimws(x)
  x1=ifelse(substr(x,1,1)=="p",word(x,2,sep="[.]+"),x)
  
  aa1=word(x1,1,sep="[0-9]+")
  aa2=word(x1,2,sep="[0-9]+")
  
  a1=convert_aa(aa1)
  a2=convert_aa(aa2)
  
  pos=gsub("[a-zA-Z*]+","",x1)
  
  paste0("p.",a1,pos,a2)
}

reverse_complement <- function(seq) {
  complement_mapping <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  seq_rev=rev(unlist(strsplit(seq,NULL)))
  paste0(unlist(lapply(seq_rev,function(x){complement_mapping[x]})),collapse = "")
}

# Convert Amino acid --------------------------------

convert_aa <- function(aa_vec) {
  # Create a lookup table of amino acid abbreviations
  aa_table <- data.frame(
    aa_3 = c("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
             "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val","*"),
    aa_1 = c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
             "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V","*"),
    stringsAsFactors = FALSE
  )
  
  # Use the lookup table to convert the 3-letter abbreviations to 1-letter abbreviations
  aa_vec_1 <- aa_table$aa_1[match(aa_vec, aa_table$aa_3)]
  
  return(aa_vec_1)
}

convert_aaPos=function(x){
  x=trimws(x)
  x1=ifelse(substr(x,1,1)=="p",word(x,2,sep="[.]+"),x)
  
  aa1=word(x1,1,sep="[0-9]+")
  aa2=word(x1,2,sep="[0-9]+")
  
  a1=convert_aa(aa1)
  a2=convert_aa(aa2)
  
  pos=gsub("[a-zA-Z*]+","",x1)
  
  paste0("p.",a1,pos,a2)
}

reverse_complement <- function(seq) {
  complement_mapping <- c("A" = "T", "T" = "A", "C" = "G", "G" = "C")
  seq_rev=rev(unlist(strsplit(seq,NULL)))
  paste0(unlist(lapply(seq_rev,function(x){complement_mapping[x]})),collapse = "")
}

