############################################# one line ############################################# 
indata=for_plot_intergrated
beta="beta"
se="se"
x="x"
x_tick_lab="x_label"
x_label="Group of discrepancy"
y_label="Disease score"

color_point="#FB2948"
color_cl="#F8AF79"
color_line="#CCA026"

forest_plot_horizontal(for_plot_disease,"beta","se","x","x_label","Group of discrepancy","Multimorbidity score")

forest_plot_horizontal=function(indata,beta,se,x,x_tick_lab,x_label,y_label){
  #dev.off()
  beta_value=as.numeric(unlist(indata[beta]))
  se_value=as.numeric(unlist(indata[se]))
  x=unlist(indata[x])
  x_tick_lab=unlist(indata[x_tick_lab])
  
  y_min=min(beta_value-se_value,na.rm = T)
  y_max=max(beta_value+se_value,na.rm = T)
  
  
  plot(NA,xlim = c(1,max(x)),ylim = c(y_min,y_max), axes = F, xlab = x_label,ylab = y_label)
  box()
  axis(side=1,at=x,labels = x_tick_lab)
  axis(side=2)
  abline(h=0,col = "lightgray", lwd = 1)
  lines(x,beta_value,col=color_line,lwd=2)
  arrows(x, beta_value+se_value, x,beta_value-se_value, angle=90, code=3, length=0.02,col=color_cl,lwd=3)
  points(x,beta_value,pch=20,col=color_point,cex=1.5)
}



############################################# multiple lines ############################################# 




######################################## Forest plot ########################################
or_male=data_or[data_or$sex=="male",] %>% mutate(color="#E4472F")
or_female=data_or[data_or$sex=="female",] %>% mutate(color="#3092CB",index=index-0.2)

data_or_new=rbind(or_male,or_female) %>% arrange(obs)

indata=data_or_new
#indata$index=nrow(indata):1

value_x_max=max(indata$ul)
value_x_min=min(indata$ll)

if(value_x_min>1){value_x_min=1}
if(value_x_max<1){value_x_max=1}

#color_in=RColorBrewer::brewer.pal(11, name = 'Spectral')[-6][1:nrow(indata)]
color_in=indata$color

pdf("2.data_temp/or.pdf",width = 5,height = 4)

par(oma=c(0,0,0,0))
par(mar=c(2,0,2,0))
par(mfcol=c(1,2)) 

plot(NA,xlim = c(-1,1),ylim = c(1,max(or_male$index)), axes = FALSE, xlab = "",ylab = "")
text(x=1,y=or_male$index,labels = or_male$Variable,cex = 0.8,pos=2)

par(mar=c(2,0.25,2,0.25))
plot(NA,xlim = c(value_x_min,value_x_max),ylim = c(1,max(indata$index)), axes = FALSE, xlab = "",ylab = "")
abline(v=1,lty=3)
points(indata$or,indata$index,pch=20,col=color_in,cex=0.5)
arrows(indata$ll, indata$index, indata$ul,indata$index, angle=90,code=3,length=0,col=color_in,lwd=2)
axis(side=1)
axis(side=2,at=or_male$index,labels=rep("",nrow(or_male)),tck=-0.025)
box()
legend("bottomright",legend=c("Male","Female"),lty=1,lwd=2,pch=20,pt.cex=0.5,col = c("#E4472F","#3092CB"),bty="n")
dev.off()
















































