label1=paste0(
  'Cutoff = ',signif(roc_result[1,8],digits = 3),'\n',
  '   SEN = ',signif(roc_result[1,6],digits = 3),'\n',
  '   SPE = ',signif(roc_result[1,7],digits = 3),'\n')

p_train=ggplot(data_train,aes(x=CT_reported_LN_status_g,y=predicted,fill=CT_reported_LN_status_g))+
  geom_dotplot(binaxis = 'y',stackdir = 'center',dotsize = 0.5)+
  geom_hline(yintercept = 0.46,colour="red")+
  xlab("CT reported LN status")+
  ylab("Predicted value")+guides(fill=F)+
  geom_text(x=2.3,y=0.5,label=label1)+
  theme_bw()+
  theme(panel.grid = element_blank())
  
tiff(paste0(outdir,"train_predicted_points.tif"))
p_train
dev.off()
