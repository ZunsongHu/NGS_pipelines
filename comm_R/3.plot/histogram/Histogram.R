draw_histogram=function(indata,var,legend_position,legend_cex,plot_title,xlabel,ymax=NA){
  value=unlist(indata[var])
  
  label1=paste0("Total N = ",length(value))
  label2=paste0("Validate N = ",length(value[!is.na(value)]))
  label3=paste0("Missing N (%) = ",length(value[is.na(value)])," (",round(length(value[is.na(value)])/length(value),4)*100,"%)")
  
  value=value[!is.na(value)]
  
  value_median=round(median(value,na.rm = T),2)
  value_q1=round(quantile(value,0.25,na.rm = T),2)
  value_q3=round(quantile(value,0.75,na.rm = T),2)
  
  legend_label=
    paste0(
      label1,"\n",
      label2,"\n",
      label3,"\n",
      
      
      paste0("Min - max: ",min(value,na.rm = T)," - ",max(value,na.rm=T)),"\n","\n",
      paste0("Median (Q1-Q3): ",value_median," (",value_q1," - ",value_q3,")"),"\n",
      paste0("Q1-1.5*IQR to Q3+1.5*IQR: ",value_q1-1.5*IQR(value),' to ',value_q3+1.5*IQR(value)," (blue)"),"\n",
      
      paste0("N outside 1.5*IQR = ",
             length(c(value[value>(value_q3+1.5*IQR(value))],
                      value[value<(value_q1-1.5*IQR(value))]))),"\n","\n",
      
      paste0("Mean (sd): ",round(mean(value,na.rm = T),2)," (",round(sd(value,na.rm = T),2),")"),"\n",
      paste0("Mean ± 3*sd: ",
             round(mean(value,na.rm = T)-3*sd(value,na.rm = T),2)," - ",
             round(mean(value,na.rm = T)+3*sd(value,na.rm = T),2)," (red)"
      ),"\n",
      paste0("N outside mean ± 3sd = ",
             length(c(value[value>(mean(value,na.rm = T)+3*sd(value,na.rm = T))],
                      value[value<(mean(value,na.rm = T)-3*sd(value,na.rm = T))])))
      
    )
  
  if(is.na(ymax)){hist(value,breaks = 100,xlab = xlabel,main=plot_title)}
  if(!is.na(ymax)){hist(value,breaks = 100,xlab = xlabel,main=plot_title,ylim=c(0,ymax))}
  
  abline(v=mean(value,na.rm = T)-3*sd(value,na.rm = T),col="red")
  abline(v=mean(value,na.rm = T)+3*sd(value,na.rm = T),col="red")
  
  abline(v=(value_q3+1.5*IQR(value)),col="blue")
  abline(v=(value_q1-1.5*IQR(value)),col="blue")
  legend(legend_position,legend=legend_label,bty="n",cex=legend_cex)
}





