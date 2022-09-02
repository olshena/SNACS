## Vanessa Kennedy single cell sequencing project

## Heatmaps: Generated under ../heatmap folder

####################################################################
####################################################################

################################################
create_heatmap=function(datObj,col_anno_var,col_dend=T,row_dend=T,hashColors=NULL,fileNameSuffix,pvSnpThres=NA) {
    
    ################################################
    ## Parameters

    subsetCellFlag="_allCells"
    row_anno_var=row_anno_name=NULL
    col_anno_name=NULL

    #pvSnpThres=NA
    sphericalK=NA
    cellClusterFileName=NA
    subsetSnpFlag=""

    cat("heatmap : pvSnpThres ",pvSnpThres,"\n")


    ################################################
    ## Parameters - optional

    dirData="../data/"

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

    ################################################
    ## Load the data
    #load(paste0(dirData,"datObj_for_heatmap.RData"))
    fName1=paste0("datObj_for_heatmap",subsetCellFlag,ifelse(is.na(sphericalK),"",paste0("_skmean",sphericalK)),ifelse(is.na(cellClusterFileName),"",paste0("_pvRanked")),ifelse(is.na(pvSnpThres),"",paste0("_pv",pvSnpThres)),"_",exptName,".RData")
    if (length(dir(dirData,pattern=fName1))==0) stop(paste0("File ",fName1," does not exist !!!"))
    load(file=paste0(dirData,fName1))

    fName=paste0(subsetCellFlag,subsetSnpFlag,fileNameSuffix)

    ################################################
    ## Folder where heatmap will be generated
    subDir="../heatmap/"
    if (!file.exists(subDir)) dir.create(file.path(subDir))

    ################################################
    ## "generate_heatmap()"  parameters
    
    if (col_dend) ncc=length(hashNames) else ncc=NA; if (!row_dend) ncr=NA else ncr=NA
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
    row_var_info=list(SnpEff_Annotation_Impact=list(color=c("brown","yellow","orange","red","white"),
                            level=c("HIGH","LOW","MODERATE","MODIFIER","")),
                            meanDepth=list(limit=c(20,80)),meanQuality=list(limit=c(25,95)))

    for (vId in which(names(col_var_info)%in%intersect(c("hash","hashCall","hashKnown","hashUnknown","hashCallMan","mclustDallCall"),names(datObj$phen)))) {
        x=sort(unique(datObj$phen[,names(col_var_info)[vId]]))
        k=match(x,col_var_info[[vId]]$level)
        k1=which(is.na(k))
        if (any(is.na(k))) {
            #col_var_info[[vId]]$level=sort(unique(datObj$phen[,names(col_var_info)[vId]]))
            #col_var_info[[vId]]$color=rainbow(length(col_var_info[[vId]]$level))
            col_var_info[[vId]]$level=c(col_var_info[[vId]]$level,x[k1])
            col_var_info[[vId]]$color=c(col_var_info[[vId]]$color,rainbow(length(k1)))
        } else {
            col_var_info[[vId]]$level=col_var_info[[vId]]$level[k]
            col_var_info[[vId]]$color=col_var_info[[vId]]$color[k]
        }
        if (!is.null(hashColors)) {
            k=match(col_var_info[[vId]]$level,names(hashColors)); k1=which(!is.na(k))
            if (length(k1)!=0) {
                k2=k[k1]
                col_var_info[[vId]]$color[k1]=hashColors[k2]
            }
        }
    }

    ## Legends
    subDir="../heatmap/legend/"
    if (!file.exists(subDir)) dir.create(file.path(subDir))
    #pdf(paste0(subDir,"cellColorLegend_hash",fileNameSuffix,".pdf"))
    pdf(paste0(subDir,"cellColorLegend_hash_",exptName,".pdf"))
    sampleColorLegend(tls=col_var_info$hashCall$level,col=col_var_info$hashCall$color)
    dev.off()
    #pdf(paste0(subDir,"cellColorLegend_mclustDallCall_",exptName,".pdf"))
    #sampleColorLegend(tls=col_var_info$mclustDallCall$level,col=col_var_info$mclustDallCall$color)
    #dev.off()
    pdf(paste0(subDir,"colorLegend_hashRaw_",exptName,".pdf"))
    heatmapColorBar(limit=col_var_info$CD45.26$limit,col=c("white","black"),main="Hash (raw)")
    dev.off()
    pdf(paste0(subDir,"cellColorLegend_cluster_skmean_",exptName,".pdf"))
    sampleColorLegend(tls=col_var_info$cluster_skmean$level,col=col_var_info$cluster_skmean$color,legendTitle="cluster_skmean")
    dev.off()
    #pdf(paste0(subDir,"snpColorLegend_SnpEff_Annotation_Impact",fileNameSuffix,".pdf"))
    pdf(paste0(subDir,"snpColorLegend_SnpEff_Annotation_Impact_",exptName,".pdf"))
    sampleColorLegend(tls=row_var_info$SnpEff_Annotation_Impact$level,col=row_var_info$SnpEff_Annotation_Impact$color)
    dev.off()
    pdf(paste0(subDir,"colorLegend_meanDepth_",exptName,".pdf"))
    heatmapColorBar(limit=row_var_info$meanDepth$limit,col=c("white","black"),main="Mean depth")
    dev.off()
    pdf(paste0(subDir,"colorLegend_meanQuality_",exptName,".pdf"))
    heatmapColorBar(limit=row_var_info$meanQuality$limit,col=c("white","black"),main="Mean quality")
    dev.off()
    pdf(paste0(subDir,"heatmap_colorLegend_",exptName,".pdf"))
    sampleColorLegend(tls=c("not mutated","mutated"),col=heatmap_color[2:1])
    dev.off()

    for (vId in which(names(col_var_info)%in%names(datObj$phen))) {
        k=match(sort(unique(datObj$phen[,names(col_var_info)[vId]])),col_var_info[[vId]]$level)
        col_var_info[[vId]]$level=col_var_info[[vId]]$level[k]
        col_var_info[[vId]]$color=col_var_info[[vId]]$color[k]
    }
    for (vId in which(names(row_var_info)%in%names(datObj$ann))) {
        k=match(sort(unique(datObj$ann[,names(row_var_info)[vId]])),row_var_info[[vId]]$level)
        row_var_info[[vId]]$level=row_var_info[[vId]]$level[k]
        row_var_info[[vId]]$color=row_var_info[[vId]]$color[k]
    }

    ################################################
    subDir="../heatmap/"
    if (subsetCellFlag!="") {subDir=paste0(subDir,subsetCellFlag,"/")} else {subDir=paste0(subDir,"_allCells/")}
    if (!file.exists(subDir)) dir.create(file.path(subDir))
    dirH=subDir

    ## Generate heatmap
    switch(outFormat,
        "png"={png(paste(dirH,"heatmap",fName,"_%1d.png",sep=""),width=480*2,height=480*2)},
        "pdf"={pdf(paste(dirH,"heatmap",fName,".pdf",sep=""),width=7,height=7)}
    )
    timeStamp=Sys.time()
    #print(format(timeStamp, "%x %X"))
    clustObj=generate_heatmap(x=datObj$mut,distfun=distfun,methodR=linkMethod,methodC=linkMethod,row_anno=hasRowAnno,col_info=datObj$phen,row_info=datObj$ann,col_dend=col_dend,row_dend=row_dend,col_lab=F,row_lab=F,zlm=limHM,ncc=ncc,ncr=ncr,heatmap_color=heatmap_color,col_anno_var=col_anno_var,col_anno_name=col_anno_name,col_var_info=col_var_info,row_anno_var=row_anno_var,row_anno_name=row_anno_name,row_var_info=row_var_info,densColor=NULL,plot_info=plotInfo,input_legend=showLegend)
    timeStamp=c(timeStamp,Sys.time())
    #print(format(timeStamp[2], "%x %X"))
    print(diff(timeStamp))
    dev.off()
    save(clustObj,file=paste0(dirH,"clustObj",fName,".RData"))
    tbl=cutCluster(clustObj$rowClust,ann=datObj$ann,nClust=ncr,rev=T)
    write.table(tbl,paste0(dirH,"clustInfoSnp",fName,".txt"), sep="\t", col.names=T, row.names=F, quote=F)
    tbl=cutCluster(clustObj$colClust,datObj$phen,nClust=ncc,rev=F)
    write.table(tbl,paste0(dirH,"clustInfoCell",fName,".txt"), sep="\t", col.names=T, row.names=F, quote=F)

}

###########################################################
###########################################################
