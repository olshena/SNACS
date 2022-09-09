createHeatmap=function(snacsObj,col_anno_var,col_anno_name=NULL,col_dend=T,row_dend=T,outputfileName="heatmap") {
    
    ################################################
    ## Parameters

    subsetCellFlag="_allCells"
    subsetCellFlag=""
    row_anno_var=row_anno_name=NULL

    cellClusterFileName=NA
    subsetSnpFlag=""

    ################################################
    ## Parameters - optional

    outFormat="png"
    outFormat="pdf"

    ## Whether to generate sample & heatmap color legends
    showLegend=T
    showLegend=F

    ################################################
    library(skmeans)
    library(RColorBrewer)
    library(marray)
    dirSrc="../src/"
    source(paste0(dirSrc,"heatmap4.R"))
    source(paste0(dirSrc,"heatmapRelated.R"))
    source(paste0(dirSrc,"generate_heatmap.R"))

    ## Distance functions available in "heatmapRelated.R"
    #getDist=getCosineDist
    #getDist=dist
    
    library(skmeans)
    getSKmeansDist=function(x) {as.dist(skmeans_xdist(x))}
    distfun=getSKmeansDist
    linkMethod="ward.D2"

    ################################################
    ## Generate heatmap of known hashes first and then all hashes
    
    if (is.null(row_anno_var)) hasRowAnno=F else hasRowAnno=T
    if (is.null(col_anno_name)) col_anno_name=col_anno_var
    if (is.null(row_anno_name)) row_anno_name=row_anno_var
    col_anno_name=paste0(col_anno_name," ") ## Cell color bar names
    row_anno_name=paste0(row_anno_name," ") ## SNP color bar names
    k=match(c("SnpEff_Annotation_Impact","meanQuality","meanDepth"),row_anno_var); k1=which(!is.na(k)); k2=k[k1]
    if (length(k1)!=0) row_anno_name[k2]=paste0(c("Impact","Quality","Depth")[k1]," ")

    #fName=paste0(subsetCellFlag,subsetSnpFlag,outputfileName)
    fName=outputfileName

    ################################################
    ## Folder where heatmap will be generated
    subDir="../output/"; if (!file.exists(subDir)) dir.create(file.path(subDir))
    subDir="../output/heatmap/"; if (!file.exists(subDir)) dir.create(file.path(subDir))

    ################################################
    ## "generate_heatmap()"  parameters
    
    if (col_dend) ncc=nrow(snacsObj$annHash) else ncc=NA; if (!row_dend) ncr=NA else ncr=NA
    heatmap_color=c("orangered1","royalblue3","white")
    limHM=c(0,1)
    #plotInfo=list(margins=c(.2,.2),cexRow=0.2,cexCol=0.2,cexRow=2,cexColSide=1.2,colorCatCol=c("skyblue", "blue", "yellow", "purple", "black", "red", "orange", "green", "cyan", "darkgreen","grey","brown","pink","salmon","palegreen"))
    #plotInfo=list(margins=c(5,.2),cexRow=0.2,cexCol=0.2,cexRow=2,cexRowSide=1,cexColSide=1.2,colorCatCol=c("skyblue", "blue", "yellow", "purple", "black", "red", "orange", "green", "cyan", "darkgreen","grey","brown","pink","salmon","palegreen"))
    plotInfo=list(margins=c(5,.2),cexRow=0.2,cexCol=0.2,cexRow=2,cexRowSide=1,cexColSide=.7,colorCatCol=c("skyblue", "blue", "yellow", "purple", "black", "red", "orange", "green", "cyan", "darkgreen","grey","brown","pink","salmon","palegreen"))
    col_var_info=list(hash=list(color=c("green3","indianred2","dodgerblue3","darkgreen","black","white"),
                            level=c("CD45-26","CD45-27","CD45-28","Multiplet","nosignal","")),
                        hashKnown=list(color=c("green3","indianred2","dodgerblue3","darkgreen","black","white"),
                            level=c("CD45-26","CD45-27","CD45-28","Multiplet","nosignal","")),
                        hashUnknown=list(color=c("green3","indianred2","dodgerblue3","darkgreen","black","white"),
                            level=c("CD45-26","CD45-27","CD45-28","Multiplet","nosignal","")),
                        hashCall=list(color=c("green3","indianred2","dodgerblue3","cyan2","darkgreen","black","white"),
                            level=c("CD45-26","CD45-27","CD45-28","Doublet","Multiplet","nosignal","")),
                        hashCallMan=list(color=c("green3","indianred2","dodgerblue3","cyan2","darkgreen","black","white"),
                            level=c("CD45-26","CD45-27","CD45-28","Doublet","Multiplet","nosignal","")),
                        mclustDallCall=list(color=c("green3","indianred2","dodgerblue3","cyan2","magenta3","yellow2","black"),
                            level=c("MixModel-26","MixModel-27","MixModel-28","MixModel-26-27","MixModel-27-28","MixModel-26-28","MixModel-nocall")),
                        CD45.26=list(limit=c(-4,4)),CD45.27=list(limit=c(-4,4)),CD45.28=list(limit=c(-4,4)),
                        cluster_hclust=list(color=c("pink","magenta","purple"),level=paste0("cluster",1:3)),
                        cluster_skmean=list(color=c("pink","magenta","purple"),level=1:3),
                        cluster_skmean2=list(color=c("pink","magenta","purple"),level=1:3),
                        cluster_skmean_skmeanSkmean=list(color=c("pink","magenta","purple"),level=1:3),
                        cluster_skmean2_skmeanSkmean=list(color=c("pink","magenta","purple"),level=1:3),
                        cluster_skmean2_hclustSkmean=list(color=c("pink","magenta","purple"),level=1:3),
                        meanDepth=list(limit=c(20,80)),meanQuality=list(limit=c(25,95)))
    col_var_info=list()
    col_var_info[["hashCall"]]=list(color=snacsObj$annHash$hashColors,level=snacsObj$annHash$hashNames)
    row_var_info=list(SnpEff_Annotation_Impact=list(color=c("brown","yellow","orange","red","white"),
                            level=c("HIGH","LOW","MODERATE","MODIFIER","")),
                            meanDepth=list(limit=c(20,80)),meanQuality=list(limit=c(25,95)))

    for (vId in which(names(col_var_info)%in%intersect(c("hash","hashCall","hashKnown","hashUnknown","hashCallMan","mclustDallCall"),names(snacsObj$annCell)))) {
        x=sort(unique(snacsObj$annCell[,names(col_var_info)[vId]]))
        k=match(x,col_var_info[[vId]]$level)
        k1=which(is.na(k))
        if (any(is.na(k))) {
            #col_var_info[[vId]]$level=sort(unique(snacsObj$annCell[,names(col_var_info)[vId]]))
            #col_var_info[[vId]]$color=rainbow(length(col_var_info[[vId]]$level))
            col_var_info[[vId]]$level=c(col_var_info[[vId]]$level,x[k1])
            col_var_info[[vId]]$color=c(col_var_info[[vId]]$color,rainbow(length(k1)))
        } else {
            col_var_info[[vId]]$level=col_var_info[[vId]]$level[k]
            col_var_info[[vId]]$color=col_var_info[[vId]]$color[k]
        }
        if (!is.null(snacsObj$annHash$hashColors)) {
            k=match(col_var_info[[vId]]$level,snacsObj$annHash$hashNames); k1=which(!is.na(k))
            if (length(k1)!=0) {
                k2=k[k1]
                col_var_info[[vId]]$color[k1]=snacsObj$annHash$hashColors[k2]
            }
        }
    }

    if (F) {
    ## Legends
    subDir="../output/heatmap/legend/"
    if (!file.exists(subDir)) dir.create(file.path(subDir))
    #pdf(paste0(subDir,"cellColorLegend_hash",outputfileName,".pdf"))
    pdf(paste0(subDir,"cellColorLegend_hash_",snacsObj$exptName,".pdf"))
    sampleColorLegend(tls=col_var_info$hashCall$level,col=col_var_info$hashCall$color)
    dev.off()
    #pdf(paste0(subDir,"cellColorLegend_mclustDallCall_",snacsObj$exptName,".pdf"))
    #sampleColorLegend(tls=col_var_info$mclustDallCall$level,col=col_var_info$mclustDallCall$color)
    #dev.off()
    pdf(paste0(subDir,"colorLegend_hashRaw_",snacsObj$exptName,".pdf"))
    heatmapColorBar(limit=col_var_info$CD45.26$limit,col=c("white","black"),main="Hash (raw)")
    dev.off()
    pdf(paste0(subDir,"cellColorLegend_cluster_skmean_",snacsObj$exptName,".pdf"))
    sampleColorLegend(tls=col_var_info$cluster_skmean$level,col=col_var_info$cluster_skmean$color,legendTitle="cluster_skmean")
    dev.off()
    #pdf(paste0(subDir,"snpColorLegend_SnpEff_Annotation_Impact",outputfileName,".pdf"))
    pdf(paste0(subDir,"snpColorLegend_SnpEff_Annotation_Impact_",snacsObj$exptName,".pdf"))
    sampleColorLegend(tls=row_var_info$SnpEff_Annotation_Impact$level,col=row_var_info$SnpEff_Annotation_Impact$color)
    dev.off()
    pdf(paste0(subDir,"colorLegend_meanDepth_",snacsObj$exptName,".pdf"))
    heatmapColorBar(limit=row_var_info$meanDepth$limit,col=c("white","black"),main="Mean depth")
    dev.off()
    pdf(paste0(subDir,"colorLegend_meanQuality_",snacsObj$exptName,".pdf"))
    heatmapColorBar(limit=row_var_info$meanQuality$limit,col=c("white","black"),main="Mean quality")
    dev.off()
    pdf(paste0(subDir,"heatmap_colorLegend_",snacsObj$exptName,".pdf"))
    sampleColorLegend(tls=c("not mutated","mutated"),col=heatmap_color[2:1])
    dev.off()
    }

    for (vId in which(names(col_var_info)%in%names(snacsObj$annCell))) {
        k=match(sort(unique(snacsObj$annCell[,names(col_var_info)[vId]])),col_var_info[[vId]]$level)
        col_var_info[[vId]]$level=col_var_info[[vId]]$level[k]
        col_var_info[[vId]]$color=col_var_info[[vId]]$color[k]
    }
    for (vId in which(names(row_var_info)%in%names(snacsObj$annSNP))) {
        k=match(sort(unique(snacsObj$annSNP[,names(row_var_info)[vId]])),row_var_info[[vId]]$level)
        row_var_info[[vId]]$level=row_var_info[[vId]]$level[k]
        row_var_info[[vId]]$color=row_var_info[[vId]]$color[k]
    }

    ################################################
    subDir="../output/heatmap/"; if (!file.exists(subDir)) dir.create(file.path(subDir))
    if (subsetCellFlag!="") {subDir=paste0(subDir,subsetCellFlag,"/"); if (!file.exists(subDir)) dir.create(file.path(subDir))}
    dirH=subDir

    ## Generate heatmap
    switch(outFormat,
        "png"={png(paste(dirH,fName,"_%1d.png",sep=""),width=480*2,height=480*2)},
        "pdf"={pdf(paste(dirH,fName,".pdf",sep=""),width=7,height=7)}
    )
    #timeStamp=Sys.time()
    ##print(format(timeStamp, "%x %X"))
    
    clustObj=generate_heatmap(x=snacsObj$mut,distfun=distfun,methodR=linkMethod,methodC=linkMethod,row_anno=hasRowAnno,col_info=snacsObj$annCell,row_info=snacsObj$annSNP,col_dend=col_dend,row_dend=row_dend,col_lab=F,row_lab=F,zlm=limHM,ncc=ncc,ncr=ncr,heatmap_color=heatmap_color,col_anno_var=col_anno_var,col_anno_name=col_anno_name,col_var_info=col_var_info,row_anno_var=row_anno_var,row_anno_name=row_anno_name,row_var_info=row_var_info,densColor=NULL,plot_info=plotInfo,input_legend=showLegend)
    #timeStamp=c(timeStamp,Sys.time())
    ##print(format(timeStamp[2], "%x %X"))
    #print(diff(timeStamp))
    dev.off()
    #save(clustObj,file=paste0(dirH,"clustObj",fName,".RData"))
    if (F) {
        tbl=cutCluster(clustObj$rowClust,ann=snacsObj$annSNP,nClust=ncr,rev=T)
        write.table(tbl,paste0(dirH,"clustInfoSnp_",fName,".txt"), sep="\t", col.names=T, row.names=F, quote=F)
        tbl=cutCluster(clustObj$colClust,ann=snacsObj$annCell,nClust=ncc,rev=F)
        write.table(tbl,paste0(dirH,"clustInfoCell_",fName,".txt"), sep="\t", col.names=T, row.names=F, quote=F)
    }
    
    invisible(clustObj)

}

###########################################################
###########################################################
