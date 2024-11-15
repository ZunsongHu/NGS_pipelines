##for Affy

liebiao<-read.table("ап╠М.txt",header=F)
for (i in 1:88){
a<-paste(liebiao[i,1])
b<-paste(liebiao[i,2])

data<-read.csv(a,as.is=T,header=T)
attach(data)
p_snp=pchisq((beta_SNP_addA1/sebeta_SNP_addA1)^2,df=1,lower=F)
p_int=pchisq((beta_SNP_bmi/sebeta_SNP_bmi)^2, df=1, lower=F)
c<-cbind(data,p_snp,p_int)
write.table(c,b,row.names=F, quote=FALSE, sep=" ") 
detach(data)
}


##Filepath:
##$ /media/000EA44100079553/xyang/GEint/addout/after_awk_out

## using R

list<-read.table("Axiom filename list.txt",header=F)

for (i in 1:116){

a<-paste(list[i,1])
b<-paste(list[i,2])


data<-read.table(a,header=T)
attach(data)
p_snp=pchisq((beta_SNP_addA1/sebeta_SNP_addA1)^2,df=1,lower=F)
p_int=pchisq((beta_SNP_bmi/sebeta_SNP_bmi)^2, df=1, lower=F)
c<-cbind(data,p_snp,p_int)
write.table(c,b,row.names=F, quote=FALSE, sep=" ") 
detach(data)
}



x<-read.table("GEint_Axiom_dbp5066_chr2_part1_add_out_new.txt",as.is=T,header=T,sep=" ")
x$beta_SNP_addA1<-as.numeric(x$beta_SNP_addA1)
x$p_snp2=pchisq((x$beta_SNP_addA1/x$sebeta_SNP_addA1)^2,df=1,lower=F)
x$p_int=pchisq((x$beta_SNP_bmi/x$sebeta_SNP_bmi)^2, df=1, lower=F)
c<-cbind(data,p_snp,p_int)



y<-read.table("GEint_map_chr1_part1.txt",as.is=T,header=T)







