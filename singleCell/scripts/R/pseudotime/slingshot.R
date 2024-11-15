# sling traj
# library(slingshot)
# 
var_cluster="group"
pseduotime_start="YN"
pseduotime_end="NY"

sce.sling=slingshot(obj_in1,clusterLabels=unlist(obj_in1[[var_cluster]]),
                    start.clus=pseduotime_start,end.clus=pseduotime_end,
                    reducedDim="PCA")

embeded=embedCurves(sce.sling,"TSNE")
embeded=slingCurves(embeded)[[1]]
embeded=data.frame(embeded$s[embeded$ord,])

sce.sling$slingPseudotime_1

plotTSNE(sce.sling,colour_by="slingPseudotime_1")+
    geom_line(data=embeded,mapping=aes(x=tSNE_1,y=tSNE_2))

