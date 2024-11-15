
########################################  Calculate annova P value of one or multiple variables ########################################
#indata=for_compare
#var=var_con_in[1]
#group="source"
#var_con_in

#cal_annova_p(snp_black,'rs1801133',"M1PFolate_result")

cal_annova_p=function(indata,group,var_list){
  
  bind_rows(lapply(var_list, function(var){
    con_group=factor(unlist(indata[group]))
    con_var=unlist(indata[var])
    data.frame(Variable=var,
               P_raw=signif(summary(aov(con_var~as.factor(con_group)))[[1]]$`Pr(>F)`[1],digits = 2),
               P=convert_p(summary(aov(con_var~as.factor(con_group)))[[1]]$`Pr(>F)`[1]),
               stringsAsFactors = F)
  }))
  
  
  
}

















