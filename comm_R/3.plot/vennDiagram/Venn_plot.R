############################################# 2 group ############################################# 
#indata=list_3006_1
#var_1="PassedQC_Previous_GWAS"
#var_2="Included_in_Dohad"
#var_3="PassedQC_Current_GWAS"

#venn_2group(list_3006_1,"PassedQC_Current_GWAS","PassedQC_Previous_GWAS",cat_name=c('Current GWAS','Previous GWAS'))

venn_2group=function(indata,var_1,var_2,cat_name,plotdir="4.out/_GWAS/10.3006list/"){
  n1=nrow(indata[unlist(indata[var_1])==1,])
  n2=nrow(indata[unlist(indata[var_2])==1,])
  n12=nrow(indata[unlist(indata[var_1])==1 & unlist(indata[var_2])==1,])
  
  grid.newpage()
  png(file=paste0(plotdir,"venn_2G_",var_1,"_",var_2,".png"),width = 600,height = 450)
  draw.pairwise.venn(n1, n2, n12, 
                     category = cat_name, 
                     fill = c("light blue", "pink"), 
                     cex=rep(2,3),
                     #alpha = rep(0.5, 2), 
                     lty = rep("blank",2),
                     cat.pos = c(0,0), 
                     cat.dist = c(0.025,0.025),
                     scaled = F,
                     cat.cex=rep(1.5,2)
                     
  )
  dev.off()
}



############################################# 3 group ############################################# 
#indata=list_3006_1
#var_1="PassedQC_Previous_GWAS"
#var_2="Included_in_Dohad"
#var_3="PassedQC_Current_GWAS"

#venn_3group(list_3006_1,"PassedQC_Previous_GWAS","Included_in_Dohad","PassedQC_Current_GWAS",
#            cat_name=c('Previous GWAS','Dohad','Current GWAS'))

venn_3group=function(indata,var_1,var_2,var_3,cat_name,plotdir="4.out/_GWAS/10.3006list/"){
  n1=nrow(indata[unlist(indata[var_1])==1,])
  n2=nrow(indata[unlist(indata[var_2])==1,])
  n3=nrow(indata[unlist(indata[var_3])==1,])
  
  n12=nrow(indata[unlist(indata[var_1])==1 & unlist(indata[var_2])==1,])
  n13=nrow(indata[unlist(indata[var_1])==1 & unlist(indata[var_3])==1,])
  n23=nrow(indata[unlist(indata[var_2])==1 & unlist(indata[var_3])==1,])
  
  n123=nrow(indata[unlist(indata[var_1])==1 & unlist(indata[var_2])==1 & unlist(indata[var_3])==1,])
  
  grid.newpage()
  png(file=paste0(plotdir,"venn_3G_",var_1,"_",var_2,"_",var_3,".png"),width = 600,height = 450)
  draw.triple.venn(area1 = n1, area2 = n2, area3 = n3, n12 = n12, n23 = n23, n13 = n13, 
                   n123 = n123,
                   category = cat_name, 
                   fill = c("skyblue", "pink1", "mediumorchid"), 
                   cex=rep(2,7),
                   #alpha = rep(0.5, 2), 
                   lty = "blank",
                   cat.dist = 0.025,
                   scaled = F,
                   cat.cex=rep(1.5,3)
  )
  dev.off()
}
