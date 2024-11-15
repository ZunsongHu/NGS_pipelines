options(warn=-1)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(ggpattern)
  library(ggrepel)
  library(ggbeeswarm)
  library(patchwork)
  library(scCustomize)
  library(ggsignif)
  library(pheatmap)
  
  
})
options(warn=0)

# remotes::install_github("trevorld/ggpattern")

#themes ----------------------
theme_paper=theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.background = element_rect(color="white",fill="white"),
                  axis.line = element_line(colour = "black"),
)

theme_boxplot=theme(panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    panel.grid.major.y = element_line(color="grey77"),
                    panel.grid.minor.y = element_line(color="grey77"),
                    panel.background = element_rect(color="white",fill="white")
)

legend_axis_text_size=function(x=14){
  theme(
    legend.text = element_text(size = x),
    axis.text =  element_text(size = x),
    axis.title = element_text(size = x)
  )
}

#rotat x axis angle --------------------
rotate_x=function(x){theme(axis.text.x = element_text(angle =x, vjust = 1, hjust=1))}


#draw_diag_barplot -----------------
draw_diag_barplot=function(diag_in,file_barplot){
  df_freq=data.frame(table(diag_in$diag)) %>% arrange(Var1) %>% arrange(Freq)
  names(df_freq)[1]="diag"
  df_freq=df_freq %>% left_join(df_ball_col) %>%
    mutate(percetage=sprintf("%.1f",100*Freq/nrow(diag_in)),
           diag_percentage=paste0(diag," (",percetage,"%)"))
  
  max_x=max(df_freq$Freq)
  
  pdf(file_barplot,width = 7,height = 14)
  par(mar=c(8,10,8,6),oma=c(6, 4.1, 4.1, 2.1))
  plot_bar=barplot(df_freq$Freq,xlab="Frequency",xaxt='n',xlim=c(0,max_x),names.arg = df_freq$Var1,col=df_freq$col,horiz=TRUE,cex.names = 1.8)
  axis(2, at=plot_bar, labels=df_freq$diag_percentage, tick=FALSE, las=2, line=-0.5, cex.axis=1.2,font=2)
  axis(1,cex.axis=1.2,font=2)
  dev.off()
}

#ggplot2 plots general --------------------------------------------------------------------------------------------------------------------------------
gg_barplot=function(x,file_out_prefix,width=5,height=10){
  x[is.na(x)]="NA";  n=length(x)
  
  df_freq=data.frame(table(x)) %>% arrange(desc(Freq)) %>%
    mutate(Group=as.character(x),
           ratio=Freq/n,
           Percentage=paste0(Group,": ",Freq," (",sprintf("%.2f",ratio*100),"% )"))
  
  ggplot(df_freq,aes(x=ratio,y=Percentage)) + geom_bar(stat="identity") +
    theme_classic()
  ggsave(paste0(file_out_prefix,".png"),width=width,height = height)
  ggsave(paste0(file_out_prefix,".pdf"),width=width,height = height)
}


gg_barplot=function(data_in,x=NULL,y=NULL,group.by=NULL,strata=NULL,cols=cols_in_default,type="identity",
                    position="dodge",bar_width=0.9,
                    y_lab=NULL,x_lab=NULL,legend_lab=NULL,
                    x_tick_label_var=NULL,
                    reverse_y=F
){
  # position in dodge,fill,stack
  col_in=cols
  
  #get lab
  if(is.null(y_lab)){y_lab=y}
  if(is.null(x_lab)){x_lab=x}
  
  if(tolower(type)=="self"){
    df_plot=data_in[c(group.by,strata)]
    names(df_plot)=c('group.by','strata')[!c(is.null(group.by),is.null(strata))]
    
    df_plot1=df_plot %>%
      filter(!is.na(group.by)) %>% mutate(n=n()) %>%
      group_by(group.by) %>% mutate(n_group.by=n()) %>%
      distinct() %>% ungroup() %>% arrange(n_group.by) %>%
      mutate(
        obs=1:n(),
        label=paste0(sprintf("%02d",obs),".",group.by," (N= ",n_group.by,")"),
        percentage=paste0(sprintf("%.1f",n_group.by/n*100),"%"),
        label_min=paste0(group.by," (",n_group.by,", ",percentage,")")
      )
  }
  
  if(tolower(type)=="identity"){
    df_plot=data_in[c(x,y,group.by,strata,x_tick_label_var)]
    names(df_plot)=c("x","y",'group.by','strata',"x_tick_label_var")[!c(is.null(x),is.null(y),is.null(group.by),is.null(strata),is.null(x_tick_label_var))]
    
    df_plot1=df_plot
  }
  
  p_base=ggplot(df_plot1,aes(x=x,y=y,fill=group.by))+
    geom_bar(stat="identity",position=position,width = bar_width)
  
  
  p1=p_base+
    scale_fill_manual(values = col_in) +
    theme(
      # legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"))+
    xlab(x_lab) +    ylab(y_lab)
  
  
  
  if(!is.null(strata)){
    p1=p1+facet_wrap(~strata,nrow = 1)
  }
  
  if(!is.null(x_tick_label_var)){
    df_x_tick_label=df_plot1 %>% ungroup() %>%  dplyr::select(x,x_tick_label_var) %>% distinct() %>% arrange(x)
    # print(df_x_tick_label)
    p1=p1+scale_x_discrete(labels=df_x_tick_label$x_tick_label_var)
  }
  
  if(!is.null(legend_lab)){p1=p1+labs(fill = legend_lab)}
  
  if(reverse_y){p1=p1 + scale_y_discrete(limits = rev(sort(unique(df_plot1$y))))}
  
  p1
  
}

gg_barplot_stack=function(df_in,var_group_x,var_group_bar,var_barLabel=NULL,var_strata=NULL,cols_fill=NULL,bar_width=0.8,
                          title=NULL,y_lab=NULL,x_lab=NULL,legend_lab=NULL,
                          x_tick_label_var=NULL,fontSize=14,x_rotate=NULL,prefix=NULL,w=5,h=6){
  
  g1=ggplot(df_in,aes(x=.data[[var_group_x]],y=per,fill=.data[[var_group_bar]]))+
    geom_bar(stat = "identity",position = "stack",width = bar_width)+
    labs(x=x_lab,y=y_lab,fill=legend_lab,title=title)+
    legend_axis_text_size(fontSize)+theme_paper
  
  if(!is.null(var_barLabel)){
    g1=g1+geom_text(aes(label=.data[[var_barLabel]]),position = position_stack(vjust = 0.5))
  }
  
  if(!is.null(cols_fill)){g1=g1+scale_fill_manual(values=cols_fill)}
  if(!is.null(x_tick_label_var)){g1=g1+scale_x_discrete(labels=x_tick_label_var)}
  if(!is.null(x_rotate)){g1=g1+rotate_x(x_rotate)}
  
  if(!is.null(var_strata)){g1=g1+facet_wrap(~.data[[var_strata]])}
  
  g1
  
  if(!is.null(prefix)){
    if(!dir.exists(dirname(prefix))){dir.create(dirname(prefix))}
    write_tsv(df_in,paste0(prefix,".df_plot.tsv"))
    ggsave(paste0(prefix,".barPlot.png"),plot = g1,width = w,height = h)
    ggsave(paste0(prefix,".barPlot.pdf"),plot = g1,width = w,height = h)
  }
  
}

draw_barplot=function(x,file_out_prefix,width=5,height=10){
  x[is.na(x)]="NA";  n=length(x)
  
  df_freq=data.frame(table(x)) %>% arrange(desc(Freq)) %>%
    mutate(Group=as.character(x),
           ratio=Freq/n,
           Percentage=paste0(Group,": ",Freq," (",sprintf("%.2f",ratio*100),"% )"))
  
  ggplot(df_freq,aes(x=ratio,y=Percentage)) + geom_bar(stat="identity") +
    theme_classic()
  ggsave(paste0(file_out_prefix,".png"),width=width,height = height)
  ggsave(paste0(file_out_prefix,".pdf"),width=width,height = height)
}

gg_boxPlot=function(df_in,var_value,var_group,var_group_fill=NULL,var_strata=NULL,size=0.6,boxColor="black",
                    highlightLevel=NULL,sizeHighlight=1.5,
                    cols=subtypeCol,
                    plot_title=NULL,x_lab=NULL,y_lab=NULL
){
  #get plot data
  df_in1=df_in[c(var_value,var_group,var_group_fill)]
  names(df_in1)=c("var_value","var_group","var_group_fill")[c(T,T,!is.null(var_group_fill))]
  
  #get labs
  if(is.null(x_lab)){x_lab=var_group}
  if(is.null(y_lab)){y_lab=var_value}
  if(is.null(plot_title)){plot_title=var_group}
  
  
  if(!is.null(var_group_fill)){
    #get color
    cols_in=cols[names(cols) %in% df_in1$var_group_fill]
    p1=ggplot(df_in1,aes(x=var_group,y=var_value,color=var_group_fill))+
      geom_boxplot(outlier.shape = NA,position=position_dodge(0.9)) +
      geom_point(position=position_dodge(1),size=size)+
      
      # geom_dotplot(binaxis='y', stackdir='centerwhole',position=position_dodge(1),dotsize=0.4)+
      # scale_color_manual(values=cols_in) +
      scale_color_manual(values=cols_in) +
      
      labs(title = plot_title, x = x_lab, y = y_lab)+
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_line(color="grey77"),
            panel.grid.minor.y = element_line(color="grey77"),
            
            panel.background = element_rect(color="white",fill="white"),
            
            axis.line = element_line(colour = "black"),
            axis.text.x = element_text(angle =45, vjust = 1, hjust=1),
            legend.position="top")
    
  } else {
    #get color
    # col_in=get_cols_cat(df_in$feature,cols_in=cols,cols_in_default)
    cols_in=cols[names(cols) %in% df_in1$var_group]
    # print(cols_in)
    
    # base boxplot Fun
    get_p=function(df){
      ggplot(df,aes(var_group,var_value))+
        geom_boxplot(outlier.shape = NA,color=boxColor) +
        geom_jitter(aes(color=var_group),size=size) +
        scale_color_manual(values=cols_in) +
        labs(title = plot_title, x = x_lab, y = y_lab)+
        theme(panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.major.y = element_line(color="grey77"),
              panel.grid.minor.y = element_line(color="grey77"),
              
              panel.background = element_rect(color="white",fill="white"),
              
              axis.line = element_line(colour = "black"),
              axis.text.x = element_text(angle =45, vjust = 1, hjust=1),
              legend.position = "none")
    }
    
    p1=get_p(df_in1)
    
    if(!is.null(highlightLevel)){
      df_highlight=df_in1[df_in1$var_group %in% highlightLevel,] %>% mutate(highlight=substr(highlightLevel,1,4))
      df_=df_in1[!df_in1$var_group %in% highlightLevel,] %>% mutate(highlight="Ref")
      
      p1=get_p(df_)+
        # geom_boxplot(data=df_highlight,aes(var_group,var_value,fill=var_group),outlier.shape = NA) +
        geom_point(data = df_highlight,aes(var_group,var_value,fill=var_group),size=sizeHighlight,color="red") +
        facet_grid(~highlight,scales = "free", space='free')
    }
  }
  p1
}


gg_boxPlot_3groups_withP=function(df_in,var_value,var_group,var_group_fill,cols_fill=NULL,x_lab=NULL,y_lab=NULL,plot_title=NULL,
                                  x_tick_label_var=NULL,fontSize=14,x_rotate=45,prefix=NULL,w=4,h=5){
  
  #get plot data
  df_in1=df_in[c(var_value,var_group,var_group_fill)]
  names(df_in1)=c("var_value","var_group","var_group_fill")[c(T,T,!is.null(var_group_fill))]
  
  #get P values
  levels_g=sort(unique(df_in1$var_group))
  if(!length(levels_g)==3){stop(paste0("Group does not include 3 levels: ", var_group))}
  
  p1=test_wilcox(df_in = df_in1[df_in1$var_group %in% c(levels_g[1:2]),],var = "var_value",g = "var_group",label = "1vs2")
  p2=test_wilcox(df_in = df_in1[df_in1$var_group %in% c(levels_g[2:3]),],var = "var_value",g = "var_group",label = "2vs3")
  p3=test_wilcox(df_in = df_in1[df_in1$var_group %in% c(levels_g[c(1,3)]),],var = "var_value",g = "var_group",label = "1vs3")
  
  df_p=bind_rows(p1,p2,p3)
  if(!is.null(prefix)){
    if(!dir.exists(dirname(prefix))){dir.create(dirname(prefix))}
    write_tsv(df_in1,paste0(prefix,".df_plot.tsv"))
    write_tsv(df_p,paste0(prefix,".df_p_wilcox.tsv"))
  }
  
  #get lab
  # if(is.null(x_lab)){x_lab=var_group}
  if(is.null(y_lab)){y_lab=var_value}
  
  #get figure
  g1=ggplot(data = df_in1,aes(x=var_group,y=var_value,fill=var_group_fill))+
    geom_boxplot(outlier.shape=NA)+
    ggbeeswarm::geom_beeswarm(data = df_in1,aes(x=var_group,y=var_value),cex = 4)+
    theme(legend.position = "none")+
    theme_paper+
    labs(title = plot_title, x = x_lab, y = y_lab)+
    legend_axis_text_size(fontSize)+
    ggsignif::geom_signif(comparisons = list(c(levels_g[1:2]),c(levels_g[2:3]),c(levels_g[c(1,3)])),
                          annotations=c(paste0("P = ",p1$P_wilcox_sci),paste0("P = ",p2$P_wilcox_sci),paste0("P = ",p3$P_wilcox_sci)),
                          step_increase = 0.12)
  
  if(!is.null(cols_fill)){g1=g1+scale_fill_manual(values=cols_fill)}
  if(!is.null(x_tick_label_var)){g1=g1+scale_x_discrete(labels=x_tick_label_var)}
  if(!is.null(x_rotate)){g1=g1+rotate_x(x_rotate)}
  
  g1
  
  if(!is.null(prefix)){
    ggsave(paste0(prefix,".boxPlot.png"),width = w,height = h)
    ggsave(paste0(prefix,".boxPlot.pdf"),width = w,height = h)
  }
}


#groups must be factor
gg_boxPlot_vlnPlot_dodge_withP=function(df_in,var_value,var_group,var_group_color,cols_color=NULL,
                                        x_lab=NULL,y_lab=NULL,plot_title=NULL,
                                        x_tick_label_var=NULL,
                                        fontSize=14,x_rotate=65,
                                        legend_pos="bottom",legend_label=NULL,legend_title=NULL,
                                        add_P=T,add_points=T,
                                        sig_bar_move=0.05,
                                        prefix=NULL,w=10,h=6
){
  
  
  
  #get plot data
  df_in1=df_in[c(var_value,var_group,var_group_color)]
  names(df_in1)=c("var_value","var_group","var_group_color")[c(T,T,!is.null(var_group_color))]
  
  if(!is.factor(df_in1$var_group)){stop("var_group is not factor!")}
  if(!is.factor(df_in1$var_group_color)){stop("var_group_color is not factor!")}
  
  #get lab
  # if(is.null(x_lab)){x_lab=var_group}
  if(is.null(y_lab)){y_lab=var_value}
  
  #get P values
  levels_g=sort(unique(df_in1$var_group))
  
  labels_=levels(df_in1$var_group)
  labels_=labels_[labels_ %in% as.character(df_in1$var_group)]
  
  df_level=data.frame(label=labels_,stringsAsFactors = F) %>% mutate(obs=1:n(),x_min=obs-0.2,x_max=obs+0.2) %>% left_join(
    df_in1 %>% group_by(var_group) %>% mutate(y_max=max(var_value),range=max(var_value)-min(var_value),ymax_1=range*sig_bar_move+y_max) %>% 
      select(var_group,y_max,range,ymax_1) %>% distinct() %>% dplyr::rename(label=var_group)
  )
  
  if(add_P){
    df_p=bind_rows(lapply(levels_g,function(level_g){
      print(level_g)
      # level_g="pDC"
      df_in2=df_in1[df_in1$var_group==level_g,]
      # print(df_in2)
      
      wilcox_out=data.frame(label = level_g,stringsAsFactors = F)
      if(nrow(df_in2)>=3 & length(unique(df_in2$var_group_color))>=2){
        wilcox_out=test_wilcox(df_in = df_in2,var = "var_value",g = "var_group_color",label = level_g)}
      
    })) %>% left_join(df_level)
    
    df_p1=df_p[df_p$P_wilcox<0.05 & !is.na(df_p$P_wilcox),]
    
    if(!is.null(prefix)){
      if(!dir.exists(dirname(prefix))){dir.create(dirname(prefix))}
      write_tsv(df_p,paste0(prefix,".df_p_wilcox.tsv"))
    }
  }
  
  if(!is.null(prefix)){
    if(!dir.exists(dirname(prefix))){dir.create(dirname(prefix))}
    write_tsv(df_in1,paste0(prefix,".df_plot.tsv"))
  }
  
  #get figure
  if(add_points){
    p1=ggplot(df_in1,aes(x=var_group,y=var_value,color=var_group_color))+
      geom_boxplot(position="dodge",outliers = F)+
      ggbeeswarm::geom_beeswarm(dodge.width=1,cex=0.8)+
      theme_paper+
      labs(title = plot_title, x = x_lab, y = y_lab,color=legend_title)+
      legend_axis_text_size(fontSize)
  }
  
  if(!add_points){
    p1=ggplot(df_in1,aes(x=var_group,y=var_value,color=var_group_color))+
      geom_boxplot(position="dodge",outliers = F)+
      theme_paper+
      labs(title = plot_title, x = x_lab, y = y_lab,color=legend_title)+
      legend_axis_text_size(fontSize)
  }
  
  if(add_P){
    if (nrow(df_p1)>=1){
      p1=p1+  geom_signif(y_position = df_p1$ymax_1, xmin = df_p1$x_min, 
                          xmax = df_p1$x_max, annotation = paste0("P = ",df_p1$P_wilcox_sci),color="black",
                          tip_length = 0)
    }
  }
  
  if(!is.null(legend_label) & !is.null(cols_color)){p1=p1+scale_color_manual(values=cols_color,labels=legend_label)} else if (!is.null(cols_color)) {p1=p1+scale_color_manual(values=cols_color)}
  if(!is.null(x_tick_label_var)){p1=p1+scale_x_discrete(labels=x_tick_label_var)}
  if(!is.null(x_rotate)){p1=p1+rotate_x(x_rotate)}
  if(!is.null(legend_pos)){p1=p1+theme(legend.position=legend_pos)}
  
  if(!is.null(prefix)){
    ggsave(paste0(prefix,".boxPlot.png"),width = w,height = h)
    ggsave(paste0(prefix,".boxPlot.pdf"),width = w,height = h)
  }
  return(p1)
}

gg_VlnPlot=function(indata,var_value,var_group,var_strata=NULL,cols_in){
  condition_=c(T,T,!is.null(var_strata))
  
  indata1=indata[c(var_value,var_group,var_strata)[condition_]]
  names(indata1)=c("var_value","var_group","var_strata")[condition_]
  
  # cols_in=cols_all[1:length(unique(indata1$var_group))]
  # names(cols_in)=unique(indata1$var_group)
  
  g1=ggplot(indata1,aes(var_value,var_group,fill=var_group))+
    geom_boxplot() +
    coord_flip() +
    scale_fill_manual(values=cols_in) +
    theme_bw()+
    theme(plot.title = element_text(face="bold", hjust = 0.5),
          axis.text.x = element_text(angle =30, vjust = 1, hjust=1),
          legend.position = "None")
  
  if(!is.null(var_strata)){g1=g1+facet_wrap(~var_strata)}
  g1
}


gg_dimPlot=function(df_in,x="uMAP_1",y="uMAP_2",var_col="diag",cols=subtypeCol,
                    size=1,axis_by=5,
                    highlightLevel=NULL,sizeHighlight=2,
                    plot_title=NULL,x_lab=NULL,y_lab=NULL,legend_title="Group",split.by=NULL){
  #get plot data
  df_in1=df_in[c(x,y,var_col,split.by)]
  if(is.null(split.by)){names(df_in1)=c("x","y","var_col")}
  if(!is.null(split.by)){names(df_in1)=c("x","y","var_col",'split.by')}
  
  #get labs
  if(is.null(x_lab)){x_lab=x}
  if(is.null(y_lab)){y_lab=y}
  if(is.null(plot_title)){plot_title=var_col}
  
  #get xaixs breaks
  min_x=min(df_in1$x);max_x=max(df_in1$x);x_breaks=round(seq(min_x,max_x,by=axis_by),0)
  min_y=min(df_in1$y);max_y=max(df_in1$y);y_breaks=round(seq(min_y,max_y,by=axis_by),0)
  
  #get color
  # col_in=get_cols_cat(df_in$feature,cols_in=cols,cols_in_default)
  col_in=cols[names(cols) %in% df_in1$var_col]
  
  p1=ggplot(df_in1,aes(x=x,y=y,col=var_col)) +
    geom_point(data = df_in1,aes(x=x,y=y,color=var_col),size=size)+
    scale_color_manual(values = col_in) +
    scale_x_continuous(breaks = x_breaks) +
    scale_y_continuous(breaks = y_breaks) +
    labs(title = plot_title, x = x_lab, y = y_lab,col=legend_title) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          axis.line = element_line(colour = "black"))+
    guides(alpha = "none")
  
  if(!is.null(highlightLevel)){
    df_highlight=df_in1[df_in1$var_col %in% highlightLevel,]
    
    p1=p1+geom_point(data=df_highlight, aes(x, y, color=var_col), size=sizeHighlight,show.legend = FALSE)+
      geom_text_repel(data=df_highlight, aes(x, y, label=var_col), fontface = "bold",show.legend = FALSE)
    
  }
  
  if(!is.null(split.by)){
    p1=p1 + facet_wrap(~split.by,scales = "free")
  }
  
  p1
}

# highlightLevel="TestSample"
gg_dimPlot_RNAseqTsne=function(df_plot,cols=subtypeCol,axis_by=10,highlightLevel=NULL,sizeHighlight=2){
  
  names(df_plot)[grepl("FeatureN",names(df_plot))]
  
  N=unique(df_plot$FeatureN)
  i=unique(df_plot$perplexityN)
  
  x_lab=paste0("tSNE dimension 1\nGene N = ", N, "; Perplexity=", i)
  y_lab="tSNE dimension 2"
  
  names(df_plot)[1:3]=c("feature","X","Y")
  
  col_in=get_cols_cat(df_plot$feature,cols_in=cols,cols_in_default)
  
  p1=ggplot() +     xlab(x_lab) +    ylab(y_lab)+    theme_bw() +
    geom_point(data=df_plot, aes(X, Y, color=feature)) +
    scale_color_manual(values = col_in) +
    scale_x_continuous(breaks = round(seq(floor(min(df_plot$X)), ceiling(max(df_plot$X)), by = axis_by),1)) +
    scale_y_continuous(breaks = round(seq(floor(min(df_plot$Y)), ceiling(max(df_plot$Y)), by = axis_by),1)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text=element_text(size=15),
          axis.title=element_text(size=15),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15,face="bold")
    )
  if(!is.null(highlightLevel)){
    df_highlight=df_plot[df_plot$feature %in% highlightLevel,]
    
    p1=p1+geom_point(data=df_highlight, aes(X, Y, color=feature), size=sizeHighlight)+
      geom_text_repel(data=df_highlight, aes(X, Y, label=feature), fontface = "bold")
    
  }
  p1
}

gg_featurePlot=function(df_in,x,y,var_col,size=1){
  df_in1=df_in[c(x,y,var_col)]
  names(df_in1)=c("x","y","var_col")
  
  ggplot(df_in1,aes(x=x,y=y,col=var_col)) +
    geom_point(size=size) +
    geom_point(data = df_in1,aes(x=x,y=y,alpha=var_col),size=size)+
    scale_colour_gradient(low='lightgrey',high='red') +
    # scale_colour_gradient2(low='lightgrey', mid = "blue",high='red') +
    labs(title = var_col, x = x, y = y,col=var_col)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    guides(alpha = "none")
  
}

# df_in=df_plot
# head(df_in)
# x="NK_percentage_reNormalization"
# y="ENSG00000164136"
# split.by="diag"
# shape.by="diag"
# size=1
#
# var="diag"
# label="group.by"
#
# extract_one_var=function(df_in,var,label){
#   if(!is.null(var)){df_out=df_in[var];names(df_out)=label}
#   if(is.null(var)){df_out=NULL}
#   df_out
# }
#
# gg_featureScatter=function(df_in,x,y,size=1,
#                            cols = "lightblue",
#                            group.by = NULL,
#                            shape.by = NULL,
#                            split.by = NULL,
#                            plot.cor = TRUE){
#
#   df_in1=df_in[c(x,y)];names(df_in1)=c('x','y')
#
#   df_in1=bind_cols(df_in1,
#                    extract_one_var(df_in,group.by,"group.by"),
#                    extract_one_var(df_in,shape.by,"shape.by"),
#                    extract_one_var(df_in,split.by,"split.by"))
#
#
#   p1=ggplot(df_in1,aes(x=x,y=y)) +
#     geom_point(col="blue",size=size) +
#     labs(x = x, y = "IL15") +
#     theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#           panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#     guides(alpha = "none")
#
#   ggsave("2042/NK_IL15.png",width=10,height=10)
#   ggsave("2042/NK_IL15.pdf",width=10,height=10)
#
#
#   if(!is.null(group.by)){p1+geom_point(data = df_in1,aes(x = x,y = y,color=group.by))}
#   if(!is.null(shape.by)){p1+geom_point(data = df_in1,aes(x = x,y = y,shape=shape.by))}
#   if(!is.null(split.by)){p2=p1+facet_wrap(~split.by)}
#
#
# ggsave("2042/split_diag_NK_IL15.png",width=10,height=10)
# ggsave("2042/split_diag_NK_IL15.pdf",width=10,height=10)




#     facet_wrap(~var_strata)+theme_bw()+
gg_geomPoint=function(df_in,x,y,size=1,x_refLine=NULL,y_refLine=NULL,x_lab=NULL,y_lab=NULL,title=NULL){
  
  df_in1=df_in[c(x,y)]
  names(df_in1)=c("x","y")[c(T,T)]
  
  if(is.null(x_lab)){x_lab=x};if(is.null(y_lab)){y_lab=y};if(is.null(title)){title=NULL};
  
  p1=ggplot(df_in1,aes(x=x,y=y)) +
    geom_point(size=size) +
    labs(x = x_lab, y = y_lab,title = title) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    guides(alpha = "none")
  
  if(is.null(y_refLine)){p1=p1}
  if(!is.null(y_refLine)){p1=p1+geom_hline(yintercept = y_refLine, linetype="dashed", color = "red")}
  
  if(is.null(x_refLine)){p1=p1}
  if(!is.null(x_refLine)){p1=p1+geom_hline(xintercept = x_refLine, linetype="dashed", color = "red")}
  
  p1
}


gg_linePlot=function(df_in,x="uMAP_1",y="uMAP_2",var_col="diag",cols=subtypeCol,
                     plot_title=NULL,x_lab=NULL,y_lab=NULL,legend_title="Group",split.by=NULL){
  
  
  
  #get plot data
  df_in1=df_in[c(x,y,var_col,split.by)]
  if(is.null(split.by)){names(df_in1)=c("x","y","var_col")}
  if(!is.null(split.by)){names(df_in1)=c("x","y","var_col",'split.by')}
  
  #get labs
  if(is.null(x_lab)){x_lab=x}
  if(is.null(y_lab)){y_lab=y}
  if(is.null(plot_title)){plot_title=var_col}
  
  
  #get color
  # col_in=get_cols_cat(df_in$feature,cols_in=cols,cols_in_default)
  col_in=cols[names(cols) %in% df_in1$var_col]
  
  p1=ggplot(df_in1,aes(x=x,y=y,col=var_col)) +
    geom_line(data = df_in1,aes(x=x,y=y,color=var_col))+
    scale_color_manual(values = col_in) +
    labs(title = plot_title, x = x_lab, y = y_lab,col=legend_title)+
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color="grey77"),
          panel.grid.minor.y = element_line(color="grey77"),
          panel.background = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          axis.line = element_line(colour = "black"))+
    guides(alpha = "none")
  
  if(!is.null(split.by)){
    p1=p1 + facet_wrap(~split.by,scales = "free")
  }
  
  p1
  
}



gg_dimHeatmap=function(indata,y_tick_lab_size=1){
  # indata=indata %>% slice_head(n=10)
  indata1=indata %>%
    mutate(id=paste0(sprintf("%08d",c(n():1)),row.names(indata))) %>% reshape2::melt( id.vars=c("id"))
  
  p1=ggplot(indata1, aes(x=variable, y=id, fill= value)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle =30, vjust = 1, hjust=1),
          axis.text.y = element_text(size=y_tick_lab_size)
    )+
    scale_y_discrete(labels=rev(row.names(indata)))
  p1
}



gg_DoHeatmap=function(){}
gg_dotPlot=function(){}
gg_elbowPlot=function(df_in,x,y,var_col,size=1){}

gg_heatmap=function(df_in,x = "labels",y = "variable",var_col = "value",
                    y_lab=NULL,x_lab=NULL,
                    y_tick_label_var=NULL,
                    reverse_y=F,
                    y_xais_side="right",
                    title="Scores (column scaled to 0-1)"
){
  
  #get data
  df_in1=df_in[c(x,y,var_col,y_tick_label_var)]
  names(df_in1)=c('x','y','var_col',"y_tick_label_var")[c(T,T,T,!is.null(y_tick_label_var))]
  
  #get lab
  if(is.null(y_lab)){y_lab=""}
  if(is.null(x_lab)){x_lab=""}
  if(is.null(title)){title=NULL}
  
  p1=ggplot(df_in1, aes(x=x, y=y, fill= var_col)) +
    geom_tile() +
    scale_fill_gradientn(colours=c("midnightblue","navyblue","deepskyblue3","lightseagreen","yellow1")) +
    xlab(x_lab) +
    ylab(y_lab) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "bottom")+
    guides(fill = guide_colourbar(ticks = FALSE,barwidth = 15,barheight = 0.8,title.position = "top",title.hjust = 0.5,title = title))
  
  if(is.null(y_tick_label_var) & reverse_y){
    p1=p1 + scale_y_discrete(position=y_xais_side,limits = rev(sort(unique(df_in1$y))))}
  
  if((!is.null(y_tick_label_var)) & reverse_y){
    df_y_tick_label=df_in1 %>% dplyr::select(y,y_tick_label_var) %>% distinct() %>% arrange(y)
    p1=p1+scale_y_discrete(position=y_xais_side,labels=rev(df_y_tick_label$y_tick_label_var),limits = rev(sort(unique(df_in1$y))))  }
  
  if((!is.null(y_tick_label_var)) & !reverse_y){
    df_y_tick_label=df_in1 %>% dplyr::select(y,y_tick_label_var) %>% distinct() %>% arrange(y)
    p1=p1+scale_y_discrete(position=y_xais_side,labels=df_y_tick_label$y_tick_label_var)
  }
  
  p1
}

gg_ridgePlot=function(df_in,x,y,fill,strata=NULL,reverse_y=F){
  df_in1=df_in[c(x,y,fill,strata)]
  names(df_in1)=c('x','y','fill',"strata")[c(TRUE,TRUE,T,!is.null(strata))]
  
  p1=ggplot(df_in1, aes(x = x, y = y, fill = fill)) +
    # geom_density_ridges( alpha = .7,scale = 0.9) +
    geom_density_ridges(
      stat = "binline", bins = 100, scale = 0.95,color=rgb(1, 1, 1, alpha=0),
      draw_baseline = FALSE
    )+
    
    scale_fill_manual(values = c("skyblue", "orange", "lightgreen")) +
    theme_ridges(grid = TRUE)+
    scale_x_continuous(breaks = c(0,0.5,1))
  
  if(!is.null(strata)){
    p1=p1+facet_wrap(~strata,nrow = 1)
  }
  
  if(reverse_y){
    p1=p1+scale_y_discrete(limits = rev(sort(unique(df_in1$y))))
  }
  
  p1
}

gg_tilePlot=function(df_in,x="obs",y="strata",var_col="group",cols=subtypeCol(),
                     y_lab=NULL,x_lab=NULL,legend_lab=NULL,reverse_y=F,
                     x_tick_label_var=NULL,
                     squared=F,add_horizontal=T,
                     add_border=F,border_lwd=1,border_linetype=1,border_color="black",
                     x.axis=T){
  #get data
  df_in1=df_in[c(x,y,var_col,x_tick_label_var)]
  names(df_in1)=c('x','y','var_col',"x_tick_label_var")[c(T,T,T,!is.null(x_tick_label_var))]
  
  #get lab
  if(is.null(y_lab)){y_lab=y}
  if(is.null(x_lab)){x_lab=x}
  
  
  
  p1=ggplot(df_in1, aes(x = x, y = y, fill = var_col)) +
    geom_tile() +
    scale_fill_manual(values = cols) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom"
    ) + xlab(x_lab) +    ylab(y_lab)
  
  if(reverse_y){p1=p1 + scale_y_discrete(limits = rev(sort(unique(df_in1$y))))}
  
  if(!is.null(legend_lab)){p1=p1+labs(fill = legend_lab)}
  
  if(!is.null(x_tick_label_var)){
    df_x_tick_label=df_in1 %>% dplyr::select(x,x_tick_label_var) %>% distinct() %>% arrange(x)
    p1=p1+scale_x_discrete(labels=df_x_tick_label$x_tick_label_var)
  }
  
  if(!x.axis){p1=p1+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  }
  
  if(add_border){p1=p1 +
    geom_tile(color = border_color,lwd = border_lwd,linetype = border_linetype)+
    theme(legend.key = element_rect(fill = NA, colour = NA))
  }
  
  if(squared){p1=p1 + coord_fixed()}
  
  if(add_horizontal){
    n1=length(unique(df_in1$y))
    n2=length(unique(df_in1$x))
    
    p1=p1+geom_line(data = data.frame(x = c(0, n2) + 0.5,
                                      y = rep(2:n1, each = 2) - 0.5),
                    aes(x = x, y = y, group = y),
                    colour="white")
  }
  p1
}









gg_circularBarplot=function(df_in,x=NULL,y,fill,N_empty_bar,strata=NULL,fill_label_var=NULL,cols=subtypeCol){
  
  #get data
  df_in1=as.data.frame(df_in[c(x,y,fill,strata,fill_label_var)])
  names(df_in1)=c('x','y','fill',"strata","fill_label_var")[c(F,T,T,!is.null(strata),!is.null(fill_label_var))]
  
  # Add lines to the initial dataset
  to_add1 = as.data.frame(matrix(NA, N_empty_bar, ncol(df_in1)))
  to_add2 = as.data.frame(matrix(NA, N_empty_bar*2, ncol(df_in1)))
  colnames(to_add1) = colnames(df_in1);colnames(to_add2) = colnames(df_in1)
  to_add1$y=0;to_add2$y=0
  df_in1_ = rbind(to_add1,df_in1,to_add2)
  df_in1_$id = seq(1, nrow(df_in1_))
  
  # get y limitation
  y_lim_upper=max(df_in1_$y,na.rm = T)
  y_lim_lower=-(y_lim_upper*0.5)
  
  #get p1
  p1=ggplot(df_in1_, aes(x = id, y = y, fill = fill)) +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values=cols)+
    ylim(y_lim_lower,y_lim_upper) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm")
    ) +
    coord_polar()
  
  
  #get df for label
  if(!is.null(fill_label_var)){
    label_data = df_in1_
    number_of_bar = nrow(label_data)
    angle = 90 - 360 * (label_data$id-1) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
    label_data$hjust = ifelse( angle < -90, 1, 0) # Get the name and the y position of each label
    label_data$angle = ifelse(angle < -90, angle+180, angle)
    
    p1=p1+  geom_text(data=label_data, aes(x=id, y=y, label=fill_label_var, hjust=hjust),angle= label_data$angle,inherit.aes = FALSE   )
    
  }
  p1
}

gg_volcanoPlot=function(df_in,P,FC,var_label,
                        P_cutoff=0.05,FC_cutoff=1,
                        x_lab=NULL,y_lab=NULL,title_=NULL,
                        col_change=c('#1C7F93', "grey77", '#782AB6'),
                        n_top_to_show=20,
                        featuresToShow=NULL,
                        point_size=1,
                        show_highlit_cicle=F,
                        show_DE_number=F,
                        point_size_highlight=3,
                        label_text_size=4,prefix=NULL,w=6,h=6,write_df_DE=T){
  
  #get label:
  if(is.null(y_lab)){y_lab=P}
  if(is.null(x_lab)){x_lab=FC}
  
  df_in1=df_in[c(P,FC,var_label)]
  names(df_in1)=c('P','FC','var_label')
  
  df_in2=df_in1 %>% tidyr::drop_na() %>% mutate(
    log_P=-log10(P),
    change=ifelse(FC< -FC_cutoff & P < P_cutoff,"downreg",
                  ifelse(FC > FC_cutoff & P < P_cutoff,"upreg","no_change")),
    group=ifelse(FC<0,"lowExp",ifelse(FC>0,"highExp",NA))
  ) %>% group_by(group) %>% arrange(P) 
  print(table(df_in2$change))
  
  if(!is.null(n_top_to_show)){
    df_in3=df_in2 %>% filter(change %in% c("downreg","upreg")) 
    
    if(nrow(df_in3)>1){
      df_in3=df_in3 %>% 
        group_by(change) %>% 
        arrange(P) %>% 
        mutate(obs=1:n(),
               label=ifelse(obs<=n_top_to_show,var_label,NA)) %>% 
        filter(!is.na(label))
    }
  }
  
  if(!is.null(featuresToShow)){
    df_in3=df_in2[unlist(df_in2['var_label']) %in% featuresToShow,] %>% 
      mutate(obs=1:n(),
             label=var_label) %>% 
      filter(!is.na(label))
  }
  
  if (!is.null(n_top_to_show) & !is.null(featuresToShow)){
    df_in3=df_in2 %>% filter(change %in% c("downreg","upreg") | var_label %in% featuresToShow)  %>% 
      group_by(change) %>% 
      arrange(P) %>% 
      mutate(obs=1:n(),
             label=ifelse(obs<=n_top_to_show,var_label,NA)) %>% 
      filter(!is.na(label))
  }
  
  # col_change=c('#1C7F93', "grey77", '#782AB6')
  if(is.null(col_change)){col_change=c('blue', "black", 'red')}
  
  names(col_change)=c("downreg","no_change","upreg")
  
  g1=ggplot(data=df_in2,aes(x=FC, y=log_P,col=change)) +
    geom_point(size=point_size) + 
    scale_color_manual(values=col_change) +
    geom_vline(xintercept=c(-FC_cutoff, FC_cutoff), col="red") +
    geom_hline(yintercept=-log10(P_cutoff), col="red")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position = "none")+
    labs(x = x_lab, y = y_lab,title = title_) 
  
  if(nrow(df_in3)>0){
    g1=g1+geom_text_repel(data = df_in3,aes(label =label), size=label_text_size, min.segment.length=0,max.overlaps = 100,colour = "black")
  }
  
  if(show_highlit_cicle){g1=g1 + geom_point(data = df_in3,size=point_size_highlight,shape =1, color="red")}
  
  if(show_DE_number){
    freqs=table(df_in2$change)
    label_change=paste0(paste0(names(freqs),": ",freqs),collapse = ", ")
    g1=g1+labs(title=paste0(title_,"\n",label_change))
  }
  
  if(!is.null(prefix)){
    if(!dir.exists(dirname(prefix))){dir.create(dirname(prefix))}
    if(write_df_DE){write_tsv(df_in2,paste0(prefix,".df_plot.tsv"))}
    ggsave(paste0(prefix,".volcanoPlot.png"),plot = g1,width = w,height = h)
    ggsave(paste0(prefix,".volcanoPlot.pdf"),plot = g1,width = w,height = h)
  }
  g1
}



gg_gseaNES=function(file_neg,file_pos,topN=20,prefix,w=20,h=10){
  df_neg=read_tsv(file_neg) %>% dplyr::arrange(desc(abs(NES))) %>% slice_head(n=topN)
  df_pos=read_tsv(file_pos) %>% dplyr::arrange(desc(abs(NES))) %>% slice_head(n=topN)
  
  starting=paste0(unique(c(word(df_neg$NAME,1,sep="_"),word(df_pos$NAME,1,sep="_"))),"_")
  
  df_gsea=bind_rows(df_neg,df_pos) %>% arrange(NES) %>% 
    mutate(
      name1=gsub(starting,"",NAME),
      name2=paste0(sprintf("%03d",n():1),"_",name1)
    )
  
  df_gsea$name2=factor(df_gsea$name2)
  
  labels_y=df_gsea$name1
  names(labels_y)=df_gsea$name2
  
  custom_theme <- theme(
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color="grey77"),
    panel.grid.minor.y = element_line(color="grey77"),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank()
  )
  
  ggplot(df_gsea,aes(x=NES,y=name2,size=SIZE))+
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey40")+
    geom_point()+
    labs(title = "GSEA analysis", x = "Normalized Enrichment Score", y = "Pathway") +
    scale_y_discrete(labels=labels_y) +
    custom_theme
  
  ggsave(paste0(prefix,".png"),width = w,height = h)
  ggsave(paste0(prefix,".pdf"),width = w,height = h)
  write_tsv(df_gsea,paste0(prefix,".tsv"))
}


gg_gseaNES_clusterProfiler=function(file_gsea=NULL,df_gsea=NULL,topN=30,prefix,w=20,h=10,title="GSEA analysis"){
  if(!is.null(file_gsea)){df_gsea=vroom::vroom(file_gsea)}
  if(is.null(df_gsea)){stop("Must provide one of gsea file or df")}
  
  df_gsea1=df_gsea  %>% 
    mutate(nes_g=ifelse(NES>0,"pos",ifelse(NES<0,"neg",NA))) %>% 
    group_by(nes_g) %>% arrange(desc(abs(NES))) %>% 
    slice_head(n=topN) %>% 
    ungroup() %>% arrange(NES) %>% 
    mutate(obs=n():1,
           name2=paste0(sprintf("%03d",n():1),"_",Description)) 
  
  df_gsea1$name2=factor(df_gsea1$name2)
  
  labels_y=df_gsea1$Description
  names(labels_y)=df_gsea1$name2
  
  custom_theme <- theme(
    panel.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color="grey77"),
    panel.grid.minor.y = element_line(color="grey77"),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, face = "italic", hjust = 0.5),
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank()
  )
  
  ggplot(df_gsea1,aes(x=NES,y=name2,size=setSize))+
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey40")+
    geom_point()+
    labs(title = title, x = "Normalized Enrichment Score", y = "Pathway") +
    scale_y_discrete(labels=labels_y) +
    custom_theme
  
  ggsave(paste0(prefix,".png"),width = w,height = h)
  ggsave(paste0(prefix,".pdf"),width = w,height = h)
  # write_tsv(df_gsea1,paste0(prefix,".tsv"))
}


# For RNAseq object -------------------------
# obj_in=obj_out
draw_DimPlot=function(obj_in,group.by="diag",cols=subtypeCol,axis_by=10,reduction="umap",
                      plot_title=NULL,
                      highlightLevel=NULL,sizeHighlight=2,size=0.8){
  
  if(!substr(reduction,1,4) %in% c("tsne","umap","pca")){stop("Reduction not defined")}
  
  df_plot=get_embeding_feature(obj_in = obj_in,feature = group.by,reduction = reduction)
  
  if(reduction=="tsne"){
    if(is.null(plot_title)){plot_title="tSNE"}
    x_lab=paste0("tSNE_1\nGene N = ", unique(df_plot$FeatureN), "; Perplexity N=", unique(df_plot$perplexityN))
    
    p1=gg_dimPlot(df_in=df_plot,x="tSNE_1",y="tSNE_2",var_col=group.by,cols=cols,
                  size=size,axis_by=axis_by,
                  highlightLevel=highlightLevel,sizeHighlight=sizeHighlight,
                  plot_title=plot_title,x_lab=x_lab,y_lab=NULL,legend_title="Group")  }
  
  if(substr(reduction,1,4)=="umap"){
    if(is.null(plot_title)){plot_title="uMAP"}
    x_lab=paste0("uMAP_1\nGene N = ", unique(df_plot$FeatureN), "; Neighbor N=", unique(df_plot$n_neighbors))
    
    p1=gg_dimPlot(df_in=df_plot,x="uMAP_1",y="uMAP_2",var_col=group.by,cols=cols,
                  size=size,axis_by=axis_by,
                  highlightLevel=highlightLevel,sizeHighlight=sizeHighlight,
                  plot_title=plot_title,x_lab=x_lab,y_lab=NULL,legend_title="Group")
  }
  p1
}



draw_featurePlot=function(obj_in,features="coverage_30X",reduction="tsne"){
  df_plot=get_embeding_feature(obj_in = obj_in,features = features,reduction = reduction)
  
  if(reduction=="tsne"){
    N=unique(df_plot$FeatureN)
    i=unique(df_plot$perplexityN)
    x_lab=paste0("tSNE dimension 1\nTop variable gene N = ", N, "; Perplexity=", i)
    y_lab="tSNE dimension 2"
    p1=gg_featurePlot(df_in = df_plot,x = 'tSNE_1',y='tSNE_2',var_col = features,size = 1) +
      xlab(x_lab) +    ylab(y_lab)
  }
  
  if(reduction=="umap"){
    N=unique(df_plot$FeatureN)
    i=unique(df_plot$n_neighbors)
    x_lab=paste0("UMAP dimension 1\nTop variable gene N = ", N, "; Neighbors=", i)
    y_lab="tSNE dimension 2"
    p1=gg_featurePlot(df_in = df_plot,x = 'uMAP_1',y='uMAP_2',var_col = features,size = 1) +
      xlab(x_lab) +    ylab(y_lab)
  }
  p1
}

feature1="ENSG00000164136"
feature2="NK_percentage_reNormalization"
shape.by="diag"
split.by="diag"

draw_featureScatter=function(obj_in,
                             feature1,
                             feature2,
                             assay_name_in = "vst",
                             group.by = NULL,
                             cols = NULL,
                             size = 1,
                             shape.by = NULL,
                             split.by=NULL,
                             plot.cor = TRUE){
  
  df_plot=get_features_df(obj_in=obj_in,assay_name_in=assay_name_in,features=unique(c(feature1,feature2,group.by,shape.by,split.by)))
  
  
}

draw_BoxPlot=function(obj_in,group.by="diag",features=c("CRLF2"),
                      cols=subtypeCol,
                      assay_name_in="vst",
                      plot_title="Box Plot",x_lab="Subtype",y_lab=NULL,
                      useGeneName=F,df_geneName=info_gtf,
                      highlightLevel=NULL,sizeHighlight=1.5){
  
  if(is.null(y_lab)){y_lab=features}
  
  if(!useGeneName){
    df_feature=get_features_df(obj_in,assay_name_in=assay_name_in,features=c(features,group.by))
    p1=gg_boxPlot(df_feature,var_value = features,var_group = group.by,plot_title = plot_title,x_lab = x_lab,y_lab = y_lab,highlightLevel = highlightLevel,cols = cols)
  }
  
  if(useGeneName){
    geneId=df_geneName$gene_id[df_geneName$gene_name %in% features][1]
    df_feature=get_features_df(obj_in,assay_name_in=assay_name_in,features=c(geneId,group.by))
    p1=gg_boxPlot(df_feature,var_value = geneId,var_group = group.by,plot_title = plot_title,x_lab = x_lab,y_lab = y_lab,highlightLevel = highlightLevel,cols = cols)
  }
  p1
}


# obj_in=obj_22_vst
# feature_panel="features_DA_postiveFC_top100"
#clustering_method %in% c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")
# annotation_names_row=T
# show_rownames=T
# var_feature_name="feature1"
colors_pheatmap <- colorRampPalette(c('darkcyan', "white", "red3"))(100)
draw_pheatmap=function(obj_in,feature_panel="features_DA",variable_n=NULL,assay_name_in="vst",
                       var_group_sample="cell_type",var_feature_name=NULL,
                       annotation_names_col=F,show_colnames=F,cluster_cols=T,
                       annotation_names_row=F,show_rownames=F,cluster_rows=T,fontsize_row=1,
                       scale="row",clustering_method="ward.D",color=colors_pheatmap){
  
  #get features name
  features_in=obj_in[[feature_panel]]
  if(!is.null(variable_n)){
    features_in=obj_in[[feature_panel]][1:min(variable_n,length((obj_in[[feature_panel]])))]
  }
  
  #get exp matirx
  df_matrix=assays(obj_in$SE[features_in,])[[assay_name_in]]
  cat("Run pheatmap: Used Feature N=", dim(df_matrix)[1], "; Used Sample N=", dim(df_matrix)[2],  "\n")
  
  #get col annotation
  df_info=as.data.frame(colData(obj_in$SE))
  df_annotation_col=df_info[var_group_sample]
  
  #get row annotation
  if(!is.null(var_feature_name)){
    df_annotation_row=rowData(obj_in$SE) %>% as.data.frame()
    df_annotation_row=df_annotation_row[features_in,]
    df_annotation_row=df_annotation_row[var_feature_name]
    names(df_annotation_row)="featureNames__"
    df_annotation_row=df_annotation_row %>% filter(!is.na(featureNames__)) %>% group_by(featureNames__) %>% mutate(n=1:n())
    feature_list=unique(df_annotation_row$featureNames__[df_annotation_row$n>=2])
    df_annotation_row=df_annotation_row %>% 
      mutate(featureNames__1=ifelse(featureNames__ %in% feature_list,paste0(featureNames__,"(",n,")"),featureNames__))
    # row.names(df_matrix)=df_annotation_row$featureNames__1
  }
  
  if (scale=="row"){
    df_matrix1=df_matrix %>% as.data.frame()
    df_matrix1$value_count=apply(df_matrix1[],1,function(x){length(unique(x))})
    df_matrix=df_matrix1[!df_matrix1$value_count==1,]
    df_matrix$value_count=NULL
    cat("After removing none variable features, becasue of row scaling: Used Feature N=", dim(df_matrix)[1], "; Used Sample N=", dim(df_matrix)[2],  "\n")
  }
  
  #draw pheatmap
  pheatmap::pheatmap(df_matrix,
                     annotation_col = df_annotation_col,color = color,
                     annotation_names_col = annotation_names_col,show_colnames = show_colnames,
                     annotation_names_row = annotation_names_row,show_rownames = show_rownames,
                     scale = scale,fontsize_row=fontsize_row,
                     cluster_rows=cluster_rows,cluster_cols = cluster_cols,
                     clustering_method=clustering_method)
  
}


draw_dimHeatmap=function(){}
draw_DoHeatmap=function(){}
draw_dotPlot=function(){}
draw_ElbowPlot=function(){}
draw_ridgePlot=function(){}
draw_VlnPlot=function(){}

draw_BarPlot=function(obj_in,group.by="diag",cols=subtypeCol){
  
  data_in=colData(obj_in$SE)[group.by]
  
  df_freq=as.data.frame(table(unlist(data_in[group.by]))) %>%
    arrange(Freq,Var1) %>%
    mutate(
      obs=1:n(),
      label=paste0(sprintf("%02d",obs),".",Var1," (N= ",Freq,")"),
      label_min=paste0(Var1," (N= ",Freq,")"))
  
  
  col_in=cols[names(cols) %in% unlist(data_in[group.by])]
  
  ggplot(df_freq,aes(x=Freq,y=label,fill=Var1))+
    geom_bar(stat="identity") +
    scale_fill_manual(values = col_in) +
    scale_y_discrete(labels= df_freq$label_min) +
    theme(legend.position = "none")
}

