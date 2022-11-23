####################################################################
####################################################################
#' Make hash calls.
#'
#' Make hash call based on mutation and hash data.
#'
#' @param snacsObj SNACSList object
#' @param prob_ecdf .
#' @param probBackgndPeak .
#' @param probBelowForeground .
#' @param propAboveBackground .
#' @param minClustSize Integer. Minimum number of cells required to be in a cluster. Default is 10
#' @param clustComparePValue P-value cutoff to compare cluster pairs
#' @return A SNACSList object
#' @export
makeHashCall=function(snacsObj,prob_ecdf=0.95,probBackgndPeak=0.6,probBelowForeground=0.75,propAboveBackground=0.5,minClustSize=10,clustComparePValue=10^-5) {

    #snacsObj; prob_ecdf=0.95; probBackgndPeak=0.6; probBelowForeground=0.75; propAboveBackground=0.5; minClustSize=10; clustComparePValue=10^-5
    
    ## -----------------------------------
    if (F) {
        subsetCellFlag=""
        
        dirData=paste0("../heatmap/",subsetCellFlag,"/")
        fName=paste0("_pvBest_ann_hashCall_mclustDallCall_",snacsObj$exptName)
        fName=paste0(subsetCellFlag,"_pvBest_ann_hashCall_",snacsObj$exptName)
        fName=paste0(subsetCellFlag,"_pvBest_ann_",snacsObj$exptName)
        load(file=paste0(dirData,"clustObj",fName,".RData"))
    }
    #clustInfo=read.table(paste0(dirData,"clustInfoCell",fName,".txt"), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=-1)
    #clustInfo=snacsObj$annCell
    clustInfo=data.frame(id=snacsObj$annCell$id,t(snacsObj$hashes),stringsAsFactors=F)
    #if (!"hashCallMan"%in%names(clustInfo)) clustInfo$hashCallMan=""
    #clustInfo$hashCall=sub("-",".",clustInfo$hashCallMan)
    clustInfo=clustInfo[match(snacsObj$hclustObj_bestSNPs$labels,clustInfo$id),]

    ## --------------------------------------
    nClust=nrow(snacsObj$annHash)
    nClustFilt=NA; clustFilt=NA; clustRename=NA
    ## --------------------------------------

    #clustInfo$clustId=clustInfo[,paste0("clustId_",nrow(snacsObj$annHash))]
    clustInfo$clustId=stats::cutree(snacsObj$hclustObj_bestSNPs,k=nrow(snacsObj$annHash))
    if (!is.na(nClustFilt[1])) {
        subClust_filt=NULL
        
        clustInfo$clustId_filt=NA
        for (k in 1:length(subClust_filt)) {
            clustInfo$clustId_filt[match(names(subClust_filt[[k]]),clustInfo$id)]=k
        }
        table(clustInfo$clustId_filt,exclude=NULL)
        clustInfo$clustId[which(clustInfo$clustId_filt%in%clustFilt)]=NA
    }
    
    hashMat=as.matrix(clustInfo[,snacsObj$annHash$hashNames])
    if (F) {
        ## Normalize hash values
        for (k in 1:ncol(hashMat)) {
            if (F) {
                xAll=hashMat[,k]; xAll[which(is.na(clustInfo$clustId))]=NA
                j=which(!is.na(clustInfo$clustId))
                x=hashMat[j,k]
                thresCenter=min(x,na.rm=T); thresScale=max(abs(x-thresCenter),na.rm=T)
                xAll[j]=((x-thresCenter)/thresScale)
                hashMat[,k]=xAll
            }
            xAll=scale(hashMat[,k])
            hashMat[,k]=xAll
        }
    }
    hashBackgnd=matrix(nrow=nrow(snacsObj$annHash),ncol=nrow(snacsObj$annHash),dimnames=list(c("mean","sd","thres"),snacsObj$annHash$hashNames))
    for (k in 1:ncol(hashMat)) {
        x=stats::density(hashMat[,k],bw="SJ",na.rm=T)
        xlim=range(x$x); ylim=range(x$y); ylim=NULL; ylim=c(0,1)
        xlim=max(abs(range(x$x))); xlim=c(-xlim,xlim)
        k1=which.max(x$y)
        if (stats::quantile(hashMat[,k],probs=probBackgndPeak)<x$x[k1]) k1=which.max(x$y[1:stats::quantile(1:(k1-1),probs=probBelowForeground)])
        xx=x$x[k1]
        x2=x$x[which(x$x<xx)]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBackgndData=x2
        hashBackgndECDF=stats::ecdf(hashBackgndData)
        #ecdfbg=1-hashBackgndECDF(x2)
        x2=hashMat[which(hashMat[,k]<xx),k]-xx; x2=c(x2,-x2); x2=x2+xx
        ecdfbg=1-hashBackgndECDF(x2)
        if (F) {
            ord=order(ecdfbg,decreasing=T); xx=x2[ord]; yy=ecdfbg[ord]
            hashBackgnd[,k]=c(mean(x2),stats::sd(x2),xx[which.min(yy>=prob_ecdf)])
        }
        hashBackgnd[,k]=c(mean(x2),stats::sd(x2),stats::quantile(x2,probs=prob_ecdf))
    }

    clustInfo$cd45Clust="Doublet"
    grp=clustInfo$clustId
    if (!is.na(clustRename[1])) {
        subClust_filt=NULL
        
        for (k in clustRename) {
            grp[match(names(subClust_filt[[k]]),clustInfo$id)]=paste0("clusterNew",k)
        }
    }
    grpUniq=sort(unique(grp))
    for (gId in 1:length(grpUniq)) {
        j=which(grp==grpUniq[gId])
        x=apply(hashMat[j,],2,mean,na.rm=T)
        clustInfo$cd45Clust[j]=names(x)[which.max(x)]
    }
    x
    
    ###########################################################
    ###########################################################
    ## Identify doublets automatically
    
    cd45Clust=clustInfo$cd45Clust
    
    clustObjThis=snacsObj$hclustObj_bestSNPs
    j=match(clustObjThis$label,clustInfo$id)
    clustInfoThis=clustInfo[j,]
    cd45MatThis=hashMat[j,]

    heightUniq=rev(unique(clustObjThis$height))
    cellId=1:nrow(clustInfoThis)
    value=rep(NA,length(cellId))
    valueOrig=rep(NA,length(cellId))
    grpHigher=stats::cutree(clustObjThis,k=1); grpHigher=grpHigher[cellId]
    for (hId in 2:(length(heightUniq)-1)) {
        atHeight=heightUniq[hId]
        hashClust=stats::cutree(clustObjThis,h=atHeight)
        grp=hashClust; grp=grp[cellId]
        
        #if (any(table(grp)<minClustSize) | all(!is.na(value))) break
        if (all(table(grp)<minClustSize) | all(!is.na(value))) break
        valueOrig=paste0(hId,"_",grp,"_orig")
        grpUniq=sort(unique(grp))
        if (length(grpUniq)>1 & any(is.na(value))) {
            for (gId1 in 1:(length(grpUniq)-1)) {
                j1=which(grp==grpUniq[gId1])
                if (length(j1)<minClustSize) next
                for (gId2 in (gId1+1):length(grpUniq)) {
                    j2=which(grp==grpUniq[gId2])
                    if (!is.na(value[j2[1]])) next
                    if (length(j2)<minClustSize) {
                        if (is.na(value[j1[1]])) value[j1]=paste0(hId,"_",gId1)
                        value[j2]=paste0(hId,"_",gId2)
                        next
                    }
                    grpThis=grp[grp%in%grpUniq[c(gId1,gId2)]]
                    j=match(names(grpThis),clustInfoThis$id)
                    #res=summary(lm(cd45MatThis[j,hashThis]~grpThis))$coef
                    res=DescTools::HotellingsT2Test(cd45MatThis[j,]~grpThis)
                    if (res$p.value>=clustComparePValue) {
                        if (is.na(value[j1[1]])) value[j1]=paste0(hId,"_",gId1)
                        value[j2]=paste0(hId,"_",gId2)
                        next
                    }
                }
            }
        }
    }
    j=which(is.na(value))
    if (length(j)!=0) value[j]=valueOrig[j]
    table(value,exclude=NULL)
    
    x=value
    x[is.na(x)]="0"
    for (k in unique(x)) {
        j=cellId[which(x==k)]
        inHashes=c()
        for (hashId in snacsObj$annHash$hashNames) {
            if (round(mean(cd45MatThis[j,hashId]>hashBackgnd["thres",hashId]),2)>=propAboveBackground) {inHashes=c(inHashes,hashId)}
        }
        if (length(inHashes)==1) {cd45Clust[j]=inHashes
        } else if (length(inHashes)==2) {cd45Clust[j]="Doublet"
        } else if (length(inHashes)>2) {cd45Clust[j]="Multiplet"
        } else {cd45Clust[j]=""}
    }

    ###########################################################
    ###########################################################
    
    snacsObj$annCell$hashCall=""
    j=match(snacsObj$hclustObj_bestSNPs$label,snacsObj$annCell$id); j1=which(!is.na(j)); j2=j[j1]
    snacsObj$annCell$hashCall[j2]=cd45Clust[j1]
    snacsObj$annCell$hashCall[snacsObj$hashCall==""]=NA
    
    cat('"hashCall" column added to "annCell" table in SNACS object\n',sep="")
    
    cat("\n")
    print(table(cd45Clust,cluster=clustInfo[match(snacsObj$hclustObj_bestSNPs$label,clustInfo$id),"clustId"],exclude=NULL,dnn=list("Hash calls",paste0("Cell clusters with best SNPs"))))

    invisible(snacsObj)
}

####################################################################
####################################################################
#' Make hash calls.
#'
#' Make hash call based on mutation and hash data.
#'
#' @param snacsObj SNACSList object
#' @param prob_ecdf .
#' @param probBackgndPeak .
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @export
generateHashDensityPlot=function(snacsObj,prob_ecdf=0.95,probBackgndPeak=0.6,outputFormat=c("","pdf","png")) {
    outputFormat=outputFormat[1]

    subsetCellFlag=""

    ## -----------------------------------
    if (outputFormat=="") {
        plotInfo=list(cexMain=1,cexLab=1,cexAxis=1)
    } else {
        dirResult="../output/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
        dirResult="../output/hashDensityPlot/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
        plotInfo=list(cexMain=2,cexLab=1.5,cexAxis=1.5)
    }
    

    ## -----------------------------------
    if (F) {
        dirData=paste0("../heatmap/",subsetCellFlag,"/")
        fName=paste0("_pvBest_ann_hashCall_mclustDallCall_",snacsObj$exptName)
        fName=paste0(subsetCellFlag,"_pvBest_ann_hashCall_",snacsObj$exptName)
        fName=paste0(subsetCellFlag,"_pvBest_ann_",snacsObj$exptName)
        load(file=paste0(dirData,"clustObj",fName,".RData"))
    }
    #clustInfo=read.table(paste0(dirData,"clustInfoCell",fName,".txt"), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=-1)
    #clustInfo=snacsObj$annCell
    clustInfo=data.frame(id=snacsObj$annCell$id,t(snacsObj$hashes),stringsAsFactors=F)
    #if (!"hashCallMan"%in%names(clustInfo)) clustInfo$hashCallMan=""
    #clustInfo$hashCall=sub("-",".",clustInfo$hashCallMan)
    clustInfo=clustInfo[match(snacsObj$hclustObj_bestSNPs$labels,clustInfo$id),]

    ## --------------------------------------
    ## For manual annotation of doublets
    nClustFilt=NA; clustFilt=NA; clustRename=NA
    ## --------------------------------------

    #clustInfo$clustId=clustInfo[,paste0("clustId_",nrow(snacsObj$annHash))]
    clustInfo$clustId=stats::cutree(snacsObj$hclustObj_bestSNPs,k=nrow(snacsObj$annHash))
    if (!is.na(nClustFilt[1])) {
        subClust_filt=NULL
        
        clustInfo$clustId_filt=NA
        for (k in 1:length(subClust_filt)) {
            clustInfo$clustId_filt[match(names(subClust_filt[[k]]),clustInfo$id)]=k
        }
        table(clustInfo$clustId_filt,exclude=NULL)
        clustInfo$clustId[which(clustInfo$clustId_filt%in%clustFilt)]=NA
    }
    
    hashMat=as.matrix(clustInfo[,snacsObj$annHash$hashNames])
    hashBackgnd=matrix(nrow=2,ncol=nrow(snacsObj$annHash),dimnames=list(c("mean","sd"),snacsObj$annHash$hashNames))
    ylimHist=rep(NA,ncol(hashMat))
    for (k in 1:ncol(hashMat)) {
        y=graphics::hist(hashMat[,k],breaks=diff(range(hashMat[,k]))/.1,plot=F)
        ylimHist[k]=max(y$counts)
    }
    ylimHist=c(0,max(ylimHist))
    for (k in 1:ncol(hashMat)) {
        switch(outputFormat,
            "png"={grDevices::png(paste0(dirResult,"densityPlot_hash",subsetCellFlag,"_",snacsObj$annHash$hashNames[k],"_",snacsObj$exptName,".png"),width=5*240,height=3*240)},
            "pdf"={grDevices::pdf(paste0(dirResult,"densityPlot_hash",subsetCellFlag,"_",snacsObj$annHash$hashNames[k],"_",snacsObj$exptName,".pdf"),width=7,height=7)}
        )
        if (outputFormat!="") graphics::par(mfcol=c(3,1))
        
        x=stats::density(hashMat[,k],bw="SJ",na.rm=T)
        xlim=range(x$x); ylim=range(x$y); ylim=NULL; ylim=c(0,1)
        xlim=max(abs(range(x$x))); xlim=c(-xlim,xlim)
        k1=which.max(x$y)
        if (stats::quantile(hashMat[,k],probs=probBackgndPeak)<x$x[k1]) k1=which.max(x$y[1:stats::quantile(1:(k1-1),probs=.75)])
        xx=x$x[k1]
        x2=x$x[which(x$x<xx)]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBGSymData=x2
        hashBackgndECDF=stats::ecdf(hashBGSymData)
        #ecdfxs=1-hashBackgndECDF(hashMat[,k])
        #ecdfbg=1-hashBackgndECDF(x2)
        x2=hashMat[which(hashMat[,k]<xx),k]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBGData=x2
        ecdfbg=1-hashBackgndECDF(hashBGData)
        vertLine=c(stats::median(x2),stats::quantile(x2,probs=prob_ecdf)); vertLineLab=c("med",prob_ecdf)
        plot(x,xlim=xlim,ylim=ylim,main=paste0(snacsObj$exptName,": ",colnames(hashMat)[k]),xlab="Hash",cex.main=plotInfo$cexMain,cex.lab=plotInfo$cexLab,cex.axis=plotInfo$cexAxis)
        graphics::lines(stats::density(x2,bw="SJ"),col="red"); graphics::abline(v=vertLine,col="green")
        plot(stats::density(x2,bw="SJ"),xlim=xlim,ylim=ylim,main="",xlab="Hash background",col="red",cex.main=plotInfo$cexMain,cex.lab=plotInfo$cexLab,cex.axis=plotInfo$cexAxis); graphics::abline(v=vertLine,col="green")
        graphics::axis(side=3,at=vertLine,labels=vertLineLab,cex.axis=plotInfo$cexAxis,las=3,col="green")
        hashBackgnd[,k]=c(mean(x2),stats::sd(x2))
        x=hashMat[,k]
        graphics::hist(x,freq=T,breaks=diff(range(x))/.1,xlim=xlim,ylim=ylimHist,main=paste0(snacsObj$exptName,": ",colnames(hashMat)[k]),xlab="Hash",ylab="Count",cex.main=plotInfo$cexMain,cex.lab=plotInfo$cexLab,cex.axis=plotInfo$cexAxis)

        if (outputFormat!="") grDevices::dev.off()
    }
}


###########################################################
###########################################################
