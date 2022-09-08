makeHashCall=function(snacsObj,prob_ecdf=0.25,probBackgndPeak=0.5,probBelowForeground=0.75,propAboveBackground=0.5,minClustSize=100,clustComparePValue=10^-5,clustAssignPValue=10^-5,hashBackgndSDThres=2) {
    
    ## -----------------------------------
    if (length(hashBackgndSDThres)==1) hashBackgndSDThres=rep(hashBackgndSDThres,nrow(snacsObj$annHash))
    hashBackgndSDThres=hashBackgndSDThres[1:nrow(snacsObj$annHash)]
    names(hashBackgndSDThres)=snacsObj$annHash$hashNames

    ## -----------------------------------
    if (F) {
        dirData=paste0("../heatmap/",subsetCellFlag,"/")
        fName=paste0("_pvBest_ann_hashCall_mclustDallCall_",snacsObj$exptName)
        fName=paste0(subsetCellFlag,"_pvBest_ann_hashCall_",snacsObj$exptName)
        fName=paste0(subsetCellFlag,"_pvBest_ann_",snacsObj$exptName)
        load(file=paste0(dirData,"clustObj",fName,".RData"))
    }
    #clustInfo=read.table(paste0(dirData,"clustInfoCell",fName,".txt"), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=-1)
    clustInfo=snacsObj$annCell
    #if (!"hashCallMan"%in%names(clustInfo)) clustInfo$hashCallMan=""
    #clustInfo$hashCall=sub("-",".",clustInfo$hashCallMan)
    clustInfo=clustInfo[match(snacsObj$hclustObj_bestSNPs$labels,clustInfo$id),]

    ## --------------------------------------
    nClust=nrow(snacsObj$annHash)
    nClustFilt=NA; clustFilt=NA; clustRename=NA
    ## --------------------------------------

    #clustInfo$clustId=clustInfo[,paste0("clustId_",nrow(snacsObj$annHash))]
    clustInfo$clustId=cutree(snacsObj$hclustObj_bestSNPs,k=nrow(snacsObj$annHash))
    if (!is.na(nClustFilt[1])) {
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
        x=density(hashMat[,k],bw="SJ",na.rm=T)
        xlim=range(x$x); ylim=range(x$y); ylim=NULL; ylim=c(0,1)
        xlim=max(abs(range(x$x))); xlim=c(-xlim,xlim)
        k1=which.max(x$y)
        if (quantile(hashMat[,k],probs=probBackgndPeak)<x$x[k1]) k1=which.max(x$y[1:quantile(1:(k1-1),probs=probBelowForeground)])
        xx=x$x[k1]
        x2=x$x[which(x$x<xx)]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBackgndData=x2
        hashBackgndECDF=ecdf(hashBackgndData)
        #ecdfbg=1-hashBackgndECDF(x2)
        x2=hashMat[which(hashMat[,k]<xx),k]-xx; x2=c(x2,-x2); x2=x2+xx
        ecdfbg=1-hashBackgndECDF(x2)
        if (F) {
            ord=order(ecdfbg,decreasing=T); xx=x2[ord]; yy=ecdfbg[ord]
            hashBackgnd[,k]=c(mean(x2),sd(x2),xx[which.min(yy>=prob_ecdf)])
        }
        hashBackgnd[,k]=c(mean(x2),sd(x2),quantile(x2,probs=prob_ecdf))
    }

    clustInfo$cd45Clust="Doublet"
    grp=clustInfo$clustId
    if (!is.na(clustRename[1])) {
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
    
    heightUniq=rev(unique(snacsObj$hclustObj_bestSNPs$height))
    atHeight=heightUniq[nrow(snacsObj$annHash)]
    hashClustMain=cutree(snacsObj$hclustObj_bestSNPs,h=atHeight)

    library(skmeans)
    getSKmeansDist=function(x) {as.dist(skmeans_xdist(x))}
    distfun=getSKmeansDist
    linkMethod="ward.D2"
    
    cd45Clust=clustInfo$cd45Clust
    for (hashThis in unique(clustInfo$cd45Clust)) {
        #hashThis=c("CD45.27")
        #hashThis=c("TS.1")
        clustObjThis=snacsObj$hclustObj_bestSNPs
        #j=match(clustObjThis$label[clustObjThis$order],clustInfo$id)
        j=match(clustObjThis$label,clustInfo$id)
        clustInfoThis=clustInfo[j,]
        cd45MatThis=hashMat[j,]
        
        heightUniq=rev(unique(clustObjThis$height))
        cellId=which(clustInfoThis$cd45Clust==hashThis)
        value=rep(NA,length(cellId))
        valueOrig=rep(NA,length(cellId))
        grpHigher=cutree(clustObjThis,k=1); grpHigher=grpHigher[cellId]
        for (hId in 2:(length(heightUniq)-1)) {
            atHeight=heightUniq[hId]
            hashClust=cutree(clustObjThis,h=atHeight)
            grp=hashClust; grp=grp[cellId]
            table(hashClust)
            x=table(grp)
            if (F) {
                if (any(x<minClustSize)) {
                    grpTmp=grp
                    for (k in which(x<minClustSize)) {
                        grpThis=grp[grp==names(x)[k]][1]
                        if (k<length(x)) {
                            grpNext=grp[grp==names(x)[k+1]][1]
                            if (grpHigher[grpThis]!=grpHigher[grpNext]) grpNext=grp[grp==names(x)[k-1]][1]
                        } else {
                            grpNext=grp[grp==names(x)[k-1]][1]
                        }
                        grpTmp[which(grp==grpThis)]=grpNext
                    }
                    grp=grpTmp
                }
                grpHigher=grp
            }
            if (any(table(grp)<minClustSize) | all(!is.na(value))) break
            if (length(x)<2) next
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
                        res=summary(lm(cd45MatThis[j,hashThis]~grpThis))$coef
                        if (res[2,4]>=clustComparePValue) {
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
            #cat("\n\n--------- ",k," ----------\n",sep="")
            j=cellId[which(x==k)]
            inHashes=c()
            for (hashId in snacsObj$annHash$hashNames) {
                #if (mean(cd45MatThis[j,hashId]>(hashBackgnd["mean",hashId]+hashBackgndSDThres[hashId]*hashBackgnd["sd",hashId]))>0.5) {inHashes=c(inHashes,hashId)}
                if (round(mean(cd45MatThis[j,hashId]>hashBackgnd["thres",hashId]),2)>=propAboveBackground) {inHashes=c(inHashes,hashId)}
            }
            if (length(inHashes)==1) {cd45Clust[j]=inHashes
            } else if (length(inHashes)==2) {cd45Clust[j]="Doublet"
            } else if (length(inHashes)>2) {cd45Clust[j]="Multiplet"
            } else {cd45Clust[j]=""}
        }
        if (F) {
            table(x)
            table(auto=cd45Clust[cellId],man=clustInfoThis$hashCall[cellId],exclude=NULL)
            table(value,man=clustInfoThis$hashCall[cellId],exclude=NULL)
            table(value,clusterOrig=clustInfoThis$cd45Clust[cellId],exclude=NULL)
            table(value,auto=cd45Clust[cellId],exclude=NULL)
        }
    }

    if (F) {
        j=match(snacsObj$hclustObj_bestSNPs$label,clustInfo$id)
        callInfo=clustInfo[j,which(!names(clustInfo)%in%colnames(hashMat))]
        nm=names(callInfo)
        callInfo=cbind(callInfo,hashMat[j,],hashCallAuto=cd45Clust,stringsAsFactors=F)
        names(callInfo)=c(nm,colnames(hashMat),"hashCallAuto")
        save(callInfo,file=paste0(dirResult,"hashCall_",snacsObj$exptName,".RData"))
    }

    ###########################################################
    ###########################################################
    ## Generate hash call table
    
    if (F) {
        dirData="../data/"
        load(paste0(dirData,"cell_barcode_",snacsObj$exptName,".RData"))
        j=match(snacsObj$hclustObj_bestSNPs$label,clustInfo$id)
        tbl=cbind(cell_barcode=clustInfo$cell_barcode[j],as.data.frame(hashMat[j,]),stringsAsFactors=F)
        j=match(tbl$cell_barcode,cell_barcode)
        callInfo=matrix(nrow=length(cell_barcode),ncol=ncol(hashMat)); colnames(callInfo)=colnames(hashMat)
        for (k in colnames(callInfo)) callInfo[j,k]=tbl[,k]
        callInfo=as.data.frame(callInfo,stringsAsFactors=F)
        callInfo=cbind(cell_barcode,callInfo,hashCall="filtered",stringsAsFactors=F)
        callInfo$hashCall[j]=cd45Clust
        callInfo$hashCall=sub(".","-",callInfo$hashCall,fixed=T)
        write.table(callInfo,paste0(dirResult,"hashCall_",snacsObj$exptName,".txt"), sep="\t", col.names=T, row.names=F, quote=F)
        rm(tbl,cell_barcode)
    }

    ###########################################################
    ###########################################################
    
    print(table(hashCall=cd45Clust,cluster=clustInfo[match(snacsObj$hclustObj_bestSNPs$label,clustInfo$id),"clustId"],exclude=NULL))
    
    snacsObj$annCell$hashCall=""
    j=match(snacsObj$hclustObj_bestSNPs$label,snacsObj$annCell$id); j1=which(!is.na(j)); j2=j[j1]
    snacsObj$annCell$hashCall[j2]=cd45Clust[j1]
    snacsObj$annCell$hashCall[snacsObj$hashCall==""]=NA

    invisible(snacsObj)
}


###########################################################
###########################################################

generateHashDensityPlot=function(snacsObj,prob_ecdf=0.25,probBackgndPeak=0.5) {
    subsetCellFlag=""
    
    ## -----------------------------------
    dirResult="../output/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
    dirResult="../output/hashDensityPlot/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))

    ## -----------------------------------
    if (F) {
        dirData=paste0("../heatmap/",subsetCellFlag,"/")
        fName=paste0("_pvBest_ann_hashCall_mclustDallCall_",snacsObj$exptName)
        fName=paste0(subsetCellFlag,"_pvBest_ann_hashCall_",snacsObj$exptName)
        fName=paste0(subsetCellFlag,"_pvBest_ann_",snacsObj$exptName)
        load(file=paste0(dirData,"clustObj",fName,".RData"))
    }
    #clustInfo=read.table(paste0(dirData,"clustInfoCell",fName,".txt"), sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=-1)
    clustInfo=snacsObj$annCell
    #if (!"hashCallMan"%in%names(clustInfo)) clustInfo$hashCallMan=""
    #clustInfo$hashCall=sub("-",".",clustInfo$hashCallMan)
    clustInfo=clustInfo[match(snacsObj$hclustObj_bestSNPs$labels,clustInfo$id),]

    ## --------------------------------------
    ## For manual annotation of doublets
    nClustFilt=NA; clustFilt=NA; clustRename=NA
    ## --------------------------------------

    #clustInfo$clustId=clustInfo[,paste0("clustId_",nrow(snacsObj$annHash))]
    clustInfo$clustId=cutree(snacsObj$hclustObj_bestSNPs,k=nrow(snacsObj$annHash))
    if (!is.na(nClustFilt[1])) {
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
        y=hist(hashMat[,k],breaks=diff(range(hashMat[,k]))/.1,plot=F)
        ylimHist[k]=max(y$counts)
    }
    ylimHist=c(0,max(ylimHist))
    for (k in 1:ncol(hashMat)) {
        png(paste0(dirResult,"densityPlot_hash",subsetCellFlag,"_",snacsObj$annHash$hashNames[k],"_",snacsObj$exptName,".png"),width=5*240,height=3*240)
        par(mfcol=c(3,1))
        x=density(hashMat[,k],bw="SJ",na.rm=T)
        xlim=range(x$x); ylim=range(x$y); ylim=NULL; ylim=c(0,1)
        xlim=max(abs(range(x$x))); xlim=c(-xlim,xlim)
        k1=which.max(x$y)
        if (quantile(hashMat[,k],probs=probBackgndPeak)<x$x[k1]) k1=which.max(x$y[1:quantile(1:(k1-1),probs=.75)])
        xx=x$x[k1]
        x2=x$x[which(x$x<xx)]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBGSymData=x2
        hashBackgndECDF=ecdf(hashBGSymData)
        #ecdfxs=1-hashBackgndECDF(hashMat[,k])
        #ecdfbg=1-hashBackgndECDF(x2)
        x2=hashMat[which(hashMat[,k]<xx),k]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBGData=x2
        ecdfbg=1-hashBackgndECDF(hashBGData)

        if (F) {
            #vertLine=c(mean(x2),mean(x2)+2*sd(x2)); vertLineLab=c("mean","2xSD")
            xx=x2; vertLine=c(mean(x2),xx[which.max(ecdfbg<0.5)],xx[which.max(ecdfbg<0.10)],xx[which.max(ecdfbg<0.05)]); vertLineLab=c("mean","50%","10%","5%")
            ord=order(ecdfbg,decreasing=T); xx=x2[ord]; yy=ecdfbg[ord]
            vertLine=c(median(xx)); vertLineLab=c("med")
            for (thres in (1-prob_ecdf)) {
                kk=c()
                k2=which(yy<=thres); if (length(k2)!=0) kk=c(kk,max(k2))
                k2=which(yy>=thres); if (length(k2)!=0) kk=c(kk,min(k2))
                if (length(kk)==2) kk=kk[1]:kk[2]
                if (length(kk)!=0) vertLine=c(vertLine,median(xx[kk])); vertLineLab=c(vertLineLab,paste0(100*thres,"%"))
            }
            vertLine=c(median(xx),xx[which.min(yy>=prob_ecdf)]); vertLineLab=c("med",prob_ecdf)
            vertLine=c(median(xx),quantile(xx,probs=prob_ecdf)); vertLineLab=c("med",prob_ecdf)
        }
        vertLine=c(median(x2),quantile(x2,probs=prob_ecdf)); vertLineLab=c("med",prob_ecdf)
        plot(x,xlim=xlim,ylim=ylim,main=paste0(snacsObj$exptName,": ",colnames(hashMat)[k]),xlab="Hash",cex.main=2,cex.lab=1.5,cex.axis=1.5)
        lines(density(x2,bw="SJ"),col="red"); abline(v=vertLine,col="green")
        plot(density(x2,bw="SJ"),xlim=xlim,ylim=ylim,main="",xlab="Hash background",col="red",cex.main=2,cex.lab=1.5,cex.axis=1.5); abline(v=vertLine,col="green")
        axis(side=3,at=vertLine,labels=vertLineLab,cex.axis=1.5,las=3,col="green")
        hashBackgnd[,k]=c(mean(x2),sd(x2))
        x=hashMat[,k]
        hist(x,freq=T,breaks=diff(range(x))/.1,xlim=xlim,ylim=ylimHist,main=paste0(snacsObj$exptName,": ",colnames(hashMat)[k]),xlab="Hash",ylab="Count",cex.main=2,cex.lab=1.5,cex.axis=1.5)
        dev.off()
    }
}


###########################################################
###########################################################
