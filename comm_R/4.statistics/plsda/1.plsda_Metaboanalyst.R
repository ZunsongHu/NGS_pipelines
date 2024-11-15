######################################## Plsda analysis using  MetaboAnalystR, with missing imputation, normaliztion and scaling ########################################
library(MetaboAnalystR)

#dir_out, directory that have the metabolites data.
#file_name, file name for metabolites data, csv type, without csv postfix.

Plsda_online=function(dir_out,file_name){
  new_dir=paste0(dir_out,file_name)
  if(!file_test("-d", new_dir)) {dir.create(new_dir)}
  
  mSet<-InitDataObjects("conc", "stat", FALSE)
  mSet<-Read.TextData(mSet, paste0(dir_out,file_name,".csv"), "rowu", "disc")
  mSet<-SanityCheckData(mSet)
  mSet<-RemoveMissingPercent(mSet, percent=0.5)
  mSet<-ImputeVar(mSet, method="knn")
  mSet<-FilterVariable(mSet, "none", "F", 25)
  mSet<-PreparePrenormData(mSet)
  
  mSet<-Normalization(mSet, "MedianNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
  mSet<-PlotNormSummary(mSet, paste0(new_dir,"/norm_0_"), "png", 72, width=NA)
  mSet<-PlotSampleNormSummary(mSet, paste0(new_dir,"/snorm_0_"), "png", 72, width=NA)
  
  
  mSet<-PLSR.Anal(mSet, reg=TRUE)
  mSet<-PLSDA.CV(mSet, "T",3, "Q2")
  mSet<-PlotPLS2DScore(mSet, paste0(new_dir,"/PLSDA_2D_Plot"), "pdf", 72, width=NA, 1,2,0.95,0,0)
  set.seed(100)
  mSet<-PLSDA.Permut(mSet, 100, "bw")
  mSet<-PlotPLS.Permutation(mSet, paste0(new_dir,"/PLSDA_permutation"), "pdf", 72, width=NA)
  #VIP
  data_vip=as.data.frame(mSet$analSet$plsda$vip.mat)
  data_vip=data_vip[1]
  names(data_vip)=paste0("VIPfor_",file_name)
  data_vip$meta_id=row.names(data_vip)
  write.csv(data_vip,paste0(new_dir,"/data_vip",".csv"),na = "",row.names = F)
  
  #Data for plsda plot
  data_for_PLSDA2DPlot=
    data.frame(
      outcome=mSet$dataSet$cls,
      score1=mSet$analSet$plsr$scores[,1],
      score2=mSet$analSet$plsr$scores[,2],
      score3=mSet$analSet$plsr$scores[,3],
      stringsAsFactors = F
    )
  write.csv(data_for_PLSDA2DPlot,paste0(new_dir,"/data_for_PLSDA2DPlot",".csv"),na = "",row.names = F)
  #dev.off()
  return(mSet)
}

######################################## plsda-plot based on metaboanalyst object ########################################


#plsda_plot_based_on_obejct(meta_obj_black,1,2,0.95,0,0)

GetColorSchema <- function(mSetObj=NA, grayscale=F){
  
  #mSetObj <- .get.mSet(mSetObj);
  lvs <- levels(mSetObj$dataSet$cls); 
  grp.num <- length(lvs);
  
  if(grayscale){
    dist.cols <- colorRampPalette(c("grey90", "grey30"))(grp.num);
  }else if(exists("colVec") && !any(colVec =="#NA")){
    dist.cols <- colVec;
  }else{
    pal18 <- c("#e6194B", "#3cb44b", "#4363d8", "#42d4f4", "#f032e6", "#ffe119", "#911eb4", "#f58231", "#bfef45",
               "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075");
    
    if(grp.num <= 18){ # update color and respect default
      dist.cols <- pal18[1:grp.num];
    }else{
      dist.cols <- colorRampPalette(pal18)(grp.num);
    }
  }
  
  colors <- vector(mode="character", length=length(mSetObj$dataSet$cls));
  for(i in 1:length(lvs)){
    colors[mSetObj$dataSet$cls == lvs[i]] <- dist.cols[i];
  }
  return (colors);
}

GetShapeSchema <- function(mSetObj=NA, show.name, grey.scale){
  if(exists("shapeVec") && all(shapeVec > 0)){
    sps <- rep(0, length=length(mSetObj$dataSet$cls));
    clsVec <- as.character(mSetObj$dataSet$cls)
    grpnms <- names(shapeVec);
    for(i in 1:length(grpnms)){
      sps[clsVec == grpnms[i]] <- shapeVec[i];
    }
    shapes <- sps;
  }else{
    if(show.name | grey.scale){
      shapes <- as.numeric(mSetObj$dataSet$cls)+1;
    }else{
      shapes <- rep(19, length(mSetObj$dataSet$cls));
    }
  }
  return(shapes);
}


plsda_plot_based_on_obejct=function (mSetObj,inx1, inx2, reg = 0.95, show = 1, grey.scale = 0, use.sparse = FALSE) {


  lv1 <- mSetObj$analSet$plsr$scores[, inx1]
  lv2 <- mSetObj$analSet$plsr$scores[, inx2]
  #xlabel <- paste("Component", inx1, "(", round(100 * mSetObj$analSet$plsr$Xvar[inx1]/mSetObj$analSet$plsr$Xtotvar, 1), "%)")
  #ylabel <- paste("Component", inx2, "(", round(100 * mSetObj$analSet$plsr$Xvar[inx2]/mSetObj$analSet$plsr$Xtotvar, 1), "%)")
  
  xlabel <- paste("Component", inx1)
  ylabel <- paste("Component", inx2)
  

  #par(mar = c(5, 5, 5, 5))
  
  text.lbls <- substr(rownames(mSetObj$dataSet$norm), 1, 12)
  
  if (mSetObj$dataSet$type.cls.lbl == "integer") {
    cls <- as.factor(as.numeric(levels(mSetObj$dataSet$cls))[mSetObj$dataSet$cls])
  }
  else {
    cls <- mSetObj$dataSet$cls
  }
  lvs <- levels(cls)
  pts.array <- array(0, dim = c(100, 2, length(lvs)))
  for (i in 1:length(lvs)) {
    inx <- mSetObj$dataSet$cls == lvs[i]
    groupVar <- var(cbind(lv1[inx], lv2[inx]), na.rm = T)
    groupMean <- cbind(mean(lv1[inx], na.rm = T), mean(lv2[inx], 
                                                       na.rm = T))
    pts.array[, , i] <- ellipse::ellipse(groupVar, centre = groupMean, 
                                         level = reg, npoints = 100)
  }
  xrg <- range(lv1, pts.array[, 1, ])
  yrg <- range(lv2, pts.array[, 2, ])
  x.ext <- (xrg[2] - xrg[1])/12
  y.ext <- (yrg[2] - yrg[1])/12
  xlims <- c(xrg[1] - x.ext, xrg[2] + x.ext)
  ylims <- c(yrg[1] - y.ext, yrg[2] + y.ext)
  cols <- GetColorSchema(mSetObj, grey.scale == 1)
  uniq.cols <- unique(cols)
  plot(lv1, lv2, xlab = xlabel, xlim = xlims, ylim = ylims, ylab = ylabel, type = "n", main = "Scores Plot")
  #grid(col = "lightgray", lty = "dotted", lwd = 1)
  grid(col = "lightgray", lty = 1, lwd = 1)
  legend.nm <- unique(as.character(sort(cls)))
  if (length(uniq.cols) > 1) {
    names(uniq.cols) <- legend.nm
  }
  for (i in 1:length(lvs)) {
    if (length(uniq.cols) > 1) {
      polygon(pts.array[, , i], col = adjustcolor(uniq.cols[lvs[i]], 
                                                  alpha = 0.2), border = NA)
    }
    else {
      polygon(pts.array[, , i], col = adjustcolor(uniq.cols, 
                                                  alpha = 0.2), border = NA)
    }
    if (grey.scale) {
      lines(pts.array[, , i], col = adjustcolor("black", 
                                                alpha = 0.5), lty = 2)
    }
  }
  pchs <- GetShapeSchema(mSetObj, show, grey.scale)
  if (grey.scale) {
    cols <- rep("black", length(cols))
  }
  if (show == 1) {
    text(lv1, lv2, label = text.lbls, pos = 4, xpd = T, 
         cex = 0.75)
    points(lv1, lv2, pch = pchs, col = cols)
  }
  else {
    if (length(uniq.cols) == 1) {
      points(lv1, lv2, pch = pchs, col = cols, cex = 1)
    }
    else {
      if (grey.scale == 1 | (exists("shapeVec") && all(shapeVec >= 
                                                       0))) {
        points(lv1, lv2, pch = pchs, col = adjustcolor(cols, 
                                                       alpha.f = 0.4), cex = 1.8)
      }
      else {
        points(lv1, lv2, pch = 21, bg = adjustcolor(cols, 
                                                    alpha.f = 0.4), cex = 2)
      }
    }
  }
  uniq.pchs <- unique(pchs)
  if (grey.scale) {
    uniq.cols <- "black"
  }
  legend("topright", legend = legend.nm, pch = uniq.pchs, bg="white",
         col = uniq.cols)
}

























































