####################################################################
####################################################################
#' Make hash calls.
#'
#' Make hash calls based on mutation and hash antibody data. The cluster information from the mutation data and the hash antibody data are used to make hash calls. Each cell is assigned one of the hash IDs or called a "Doublet" when it is a mixture of hashes. The distribution of the antibody measure of a hash is expected to be a mixture of a background and a foreground signal where the latter signals the presence of antibody for the hash. Expecting the mutation data of cells with similar hash profiles to cluster together, the cells are split into sub-clusters and each cluster is assigned a hash ID if enough proportion of cells in the cluster have signal above the background for that hash. The empirical background distribution of a hash antibody is estimated by taking the background mode of the data and generating an empirical symmetric distribution around it. The sub-clusters are created using the cell clusters from the second round of hierarchical clustering of the mutation data. First the antibody measures of the top two cluster pairs are compared using Hotelling's T2 test. If there is significant difference , the two clusters are splits into sub-clusters and each of those cluster pairs are tested for difference. The tree is traversed until there are no pairs of significantly different cluster pairs. The splitting of a cluster stops it reaches a minimum number of cells.
#'
#' @param snacsObj SNACSList object
#' @param backgndThreshold Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell. Default is 0.95. Range is 0-1
#' @param cellProportionAboveBackgnd Proportion of cells in a cluster that have to have signal above the background of a hash, determined by "backgndThreshold" parameter, for the cluster to be assigned to the hash. Default is 0.5. Range is 0-1
#' @param cellProportionBelowBackgndMode Maximum proportion of cells which can be below the mode of the estimated hash background distribution. Default is 0.6. Range is 0-1
#' @param cellProportionForModeDetection Proportion of cells used to estimate mode of the background distribution. Used only if "cellProportionBelowBackgndMode" threshold is not met; otherwise, all cells are used. Default is 0.75. Range is 0-1
#' @param minClustSize Minimum number of cells required to be in a cluster. Default is 2
#' @param minClustSizeRnd2 Minimum number of cells required to be in a cluster for making second round of hash calls. Inf means the second round is not run. Default is Inf
#' @param backgndThresRnd2 Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell for making second round of hash calls. Default is 0.75. Range is 0-1
#' @param clustComparePValue P-value threshold to compare cluster pairs
#' @param maxClustSampleSize Maximum number of cells in a cluster that will be used to compare. If more, then this number of cells will be sampled. Default is 5000
#' @param clustCompareMethod Test used to compare clusters. Default is t-test
#' @return A SNACSList object
#' @export
makeHashCall=function(snacsObj,backgndThreshold=0.95,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,cellProportionAboveBackgnd=0.5,minClustSize=2,minClustSizeRnd2=Inf,backgndThresRnd2=0.75,clustComparePValue=10^-5,maxClustSampleSize=Inf,clustCompareMethod=c("t","hotelling")) {

    ## -----------------------------------
    clustCompareMethod=clustCompareMethod[1]
    ## -----------------------------------
    
    clustInfo=data.frame(id=snacsObj$annCell$id,t(snacsObj$hashes),stringsAsFactors=F)
    clustInfo=clustInfo[match(snacsObj$hclustObj_bestSNPs$labels,clustInfo$id),]

    ## --------------------------------------
    ## For manual annotation of doublets
    nClust=nrow(snacsObj$annHash)
    nClustFilt=NA; clustFilt=NA; clustRename=NA
    ## --------------------------------------

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
    hashBackgnd=matrix(nrow=4,ncol=nrow(snacsObj$annHash),dimnames=list(c("mean","sd","thres","thresRnd2"),snacsObj$annHash$hashNames))
    for (k in 1:ncol(hashMat)) {
        x=stats::density(hashMat[,k],bw="SJ",na.rm=T)
        xlim=range(x$x); ylim=range(x$y); ylim=NULL; ylim=c(0,1)
        xlim=max(abs(range(x$x))); xlim=c(-xlim,xlim)
        k1=which.max(x$y)
        if (stats::quantile(hashMat[,k],probs=cellProportionBelowBackgndMode)<x$x[k1]) k1=which.max(x$y[1:stats::quantile(1:(k1-1),probs=cellProportionForModeDetection)])
        xx=x$x[k1]
        x2=x$x[which(x$x<xx)]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBackgndData=x2
        hashBackgndECDF=stats::ecdf(hashBackgndData)
        x2=hashMat[which(hashMat[,k]<xx),k]-xx; x2=c(x2,-x2); x2=x2+xx
        ecdfbg=1-hashBackgndECDF(x2)
        hashBackgnd[,k]=c(mean(x2),stats::sd(x2),stats::quantile(x2,probs=c(backgndThreshold,backgndThresRnd2)))
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
                    if (length(j)>maxClustSampleSize) {
                        set.seed(12345); j=sample(1:length(grpThis),maxClustSampleSize,replace=FALSE)
                        grpThis=grpThis[j]
                    }
                    j=match(names(grpThis),clustInfoThis$id)
                    if (clustCompareMethod=="hotelling") {
                        res=DescTools::HotellingsT2Test(cd45MatThis[j,]~grpThis)
                        pv=res$p.value
                    } else {
                        pv=rep(NA,ncol(cd45MatThis))
                        for (k in 1:ncol(cd45MatThis)) {
                            res=stats::t.test(cd45MatThis[j,k]~grpThis,var.equal=T)
                            pv[k]=res$p.value
                        }
                        pv=min(pv*ncol(cd45MatThis))
                    }
                    if (pv>=clustComparePValue) {
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
            if (round(mean(cd45MatThis[j,hashId]>hashBackgnd["thres",hashId]),2)>=cellProportionAboveBackgnd) {inHashes=c(inHashes,hashId)}
        }
        if (F) {
            if (length(inHashes)==1) {cd45Clust[j]=inHashes
            } else if (length(inHashes)==2) {cd45Clust[j]="Doublet"
            } else if (length(inHashes)>2) {cd45Clust[j]="Multiplet"
            } else {cd45Clust[j]=""}
        }
        cd45Clust[j]=paste(inHashes,collapse="_")
    }
    
    ## Mark cluster splits
    j=clustObjThis$order
    grp=as.integer(as.factor(x[j]))
    j=which(diff(grp)!=0)
    jj=j
    j1=c(1,j+1); j2=c(j,length(grp))
    for (j in 1:length(j1)) {
        if (any(grp[j1[j]:j2[j]]!=grp[j1[j]])) print(j)
    }
    clustSplit=rep(1,length(grp))
    for (j in 1:length(j1)) {
        if (j%%2==0) clustSplit[j1[j]:j2[j]]=2
    }
    clustSplit=clustSplit[match(clustInfo$id,clustObjThis$label[clustObjThis$order])]
    subClust=x
    
    ## Distance to centroid
    grp=cd45Clust
    grpUniq=sort(unique(grp))
    centroid=matrix(nrow=ncol(cd45MatThis),ncol=length(grpUniq),dimnames=list(colnames(cd45MatThis),grpUniq))
    for (gId in 1:length(grpUniq)) {
        j=which(grp==grpUniq[gId])
        centroid[,gId]=apply(cd45MatThis[j,],2,mean,na.rm=T)
    }
    dist2centroidMat=matrix(nrow=length(grpUniq),ncol=nrow(cd45MatThis),dimnames=list(grpUniq,rownames(cd45MatThis)))
    j=which(grp%in%colnames(cd45MatThis))
    j=1:nrow(cd45MatThis)
    for (gId in 1:length(grpUniq)) {
        dist2centroidMat[gId,j]=apply(cd45MatThis[j,],1,function(x,centroid) {sum((x-centroid)^2)},centroid=centroid[,gId])
    }
    if (F) {
        ## Singleton based on distance to centroid
        
        singletonThreshold=Inf

        cd45Clust2=cd45Clust
        if (all(singletonThreshold>=0) & all(singletonThreshold<=1)) {
            thres=rep(NA,length(snacsObj$annHash$hashNames)); names(thres)=snacsObj$annHash$hashNames
            for (hashId in snacsObj$annHash$hashNames) {
                j=which(grp==hashId)
                thres[hashId]=stats::quantile(dist2centroidMat[hashId,j],probs=singletonThreshold[2])
            }
            hashId2=grpUniq[!grpUniq%in%snacsObj$annHash$hashNames]
            thresDoub=rep(NA,length(hashId2)); names(thresDoub)=hashId2
            for (hashId in names(thresDoub)) {
                j=which(grp==hashId)
                thresDoub[hashId]=stats::quantile(dist2centroidMat[hashId,j],probs=singletonThreshold[1])
            }
            grp2=cd45Clust
            for (hashId in snacsObj$annHash$hashNames) {
                if (F) {
                    j=which(grp==hashId)
                    for (hashId2 in names(thresDoub)) {
                        thresDoub[hashId]=stats::quantile(dist2centroidMat[hashId2,j],probs=singletonThreshold[1])
                    }
                }
                j1=which(dist2centroidMat[hashId,]>thres[hashId] & grp==hashId)
                #j1=which(grp==hashId)
                if (length(j1)!=0) {
                    inHashes=c()
                    for (hashId2 in names(thresDoub)) {
                        j=which(dist2centroidMat[hashId2,j1]<=thresDoub[hashId2])
                        if (length(j)!=0) {
                            j=j1[j]
                            inHashes=c(inHashes,strsplit(hashId2,"_")[[1]])
                        }
                        grp2[j]=paste(sort(unique(inHashes)),collapse="_")
                    }
                }
            }
            cd45Clust2=grp2
            #table(cd45Clust,cd45Clust2,exclude=NULL)
        }
    }
    
    if (F) {
        ## Singleton based on hash value
        if (all(singletonThreshold>=0) & all(singletonThreshold<=1)) {
            thres=matrix(rep(NA,2*nrow(snacsObj$annHash)),ncol=2); rownames(thres)=snacsObj$annHash$hashNames; colnames(thres)=c("high","low")
            for (hashId in snacsObj$annHash$hashNames) {
                j=which(cd45MatThis[,hashId]>hashBackgnd["thres",hashId])
                thres[hashId,"high"]=stats::quantile(cd45MatThis[j,hashId],probs=singletonThreshold["high"])
                thres[hashId,"low"]=stats::quantile(cd45MatThis[j,hashId],probs=singletonThreshold["low"])
            }
            for (hashId in snacsObj$annHash$hashNames) {
                j=which(cd45MatThis[,hashId]>=thres[hashId,"high"])
                for (hashId2 in snacsObj$annHash$hashNames[snacsObj$annHash$hashNames!=hashId]) {
                    j=j[which(cd45MatThis[j,hashId2]<thres[hashId2,"low"])]
                }
                if (length(j)!=0) {
                    #for (jj in j) cd45Clust[jj]=paste(sort(unique(c(strsplit(cd45Clust[jj],"_")[[1]],hashId))),collapse="_")
                    cd45Clust[j]=hashId
                }
            }
        }
    }

    ## Narrow regions
    if (is.finite(minClustSizeRnd2)) {
        cellId=which(cd45Clust%in%snacsObj$annHash$hashNames)
        if (F) {
            ## Splits based on hclust

            value=rep(NA,length(cellId))
            grpHigher=stats::cutree(clustObjThis,k=1); grpHigher=grpHigher[cellId]
            for (hId in 2:(length(heightUniq)-1)) {
                atHeight=heightUniq[hId]
                hashClust=stats::cutree(clustObjThis,h=atHeight)
                grp=hashClust; grp=grp[cellId]
                x=table(grp)
                if (all(x>minClustSizeRnd2)) next
                if (all(x<minClustSize)) break
                grpUniq=sort(unique(grp))
                #if (length(grpUniq)>1 & any(is.na(value))) {
                if (length(grpUniq)>1) {
                    for (gId1 in 1:(length(grpUniq)-1)) {
                        j1=which(grp==grpUniq[gId1])
                        if (length(j1)<minClustSize) next
                        #if (length(j1)<minClustSizeRnd2) next
                        for (gId2 in (gId1+1):length(grpUniq)) {
                            j2=which(grp==grpUniq[gId2])
                            #if (!is.na(value[j2[1]])) next
                            if (length(j2)<minClustSize) {
                            #if (length(j2)<minClustSizeRnd2) {
                                #if (is.na(value[j1[1]])) value[j1]=paste0(hId,"_",gId1)
                                #value[j2]=paste0(hId,"_",gId2)
                                next
                            }
                            grpThis=grp[grp%in%grpUniq[c(gId1,gId2)]]
                            if (length(j)>maxClustSampleSize) {
                                set.seed(12345); j=sample(1:length(grpThis),maxClustSampleSize,replace=FALSE)
                                grpThis=grpThis[j]
                            }
                            j=match(names(grpThis),clustInfoThis$id)
                            if (clustCompareMethod=="hotelling") {
                                res=DescTools::HotellingsT2Test(cd45MatThis[j,]~grpThis)
                                pv=res$p.value
                            } else {
                                pv=rep(NA,ncol(cd45MatThis))
                                for (k in 1:ncol(cd45MatThis)) {
                                    res=stats::t.test(cd45MatThis[j,k]~grpThis,var.equal=T)
                                    pv[k]=res$p.value
                                }
                                pv=min(pv*ncol(cd45MatThis))
                            }
                            #if (pv>=clustComparePValue) {
                            if (pv<clustComparePValue) {
                                #if (is.na(value[j1[1]])) value[j1]=paste0(hId,"_",gId1)
                                value[j1]=paste0(hId,"_",gId1)
                                value[j2]=paste0(hId,"_",gId2)
                                next
                            }
                        }
                    }
                }
            }
            #valueOrig=paste0(hId,"_",grp,"_orig")
            #j=which(is.na(value))
            #if (length(j)!=0) value[j]=valueOrig[j]
            table(value,exclude=NULL)
            x=value
        }

        if (F) {
            ## Splits based on genotype
            x=as.character(snacsObj$mut[1,])
            if (nrow(snacsObj$mut)>1) {
                for (i in 2:nrow(snacsObj$mut)) {
                    x=paste(x,snacsObj$mut[i,])
                }
            }
            x=x[cellId]
        }
        if (T) {
            ## Splits based on cbs
            #x1=cd45MatThis[clustObjThis$order,]
            x1=t(dist2centroidMat[,clustObjThis$order])
            x2 <- DNAcopy::CNA(genomdat=x1,
            chrom=rep(1,nrow(x1)),maploc=1:nrow(x1),
                              data.type="logratio",sampleid=colnames(x1))
            x3 <- DNAcopy::segment(x2,alpha=0.1,verbose=0)
            x3=x3$output
            x=rep(NA,nrow(x1))
            names(x)=clustObjThis$label[clustObjThis$order]
            grp=cd45Clust[clustObjThis$order]
            grpUniq=snacsObj$annHash$hashNames
            for (gId in 1:length(grpUniq)) {
                j1=range(which(grp==grpUniq[gId]))
                x2=rep(NA,nrow(x1))
                p1=1
                for (k in which(x3$ID==grpUniq[gId])) {
                    x2[x3$loc.start[k]:x3$loc.end[k]]=p1
                    p1=p1+1
                }
                x[j1[1]:j1[2]]=paste0(grpUniq[gId],"_",x2[j1[1]:j1[2]])
            }
            x=x[match(clustObjThis$label,names(x))]
            x=x[cellId]

        }

        grp=cd45Clust
        #x[is.na(x)]="0"
        for (k in unique(x[!is.na(x)])) {
            j=cellId[which(x==k)]
            if (length(j)>minClustSizeRnd2) next
            inHashes=c()
            for (hashId in snacsObj$annHash$hashNames) {
                if (round(mean(cd45MatThis[j,hashId]>hashBackgnd["thresRnd2",hashId]),2)>=cellProportionAboveBackgnd) {inHashes=c(inHashes,hashId)}
            }
            #if (length(inHashes)>1) grp[j]=paste(inHashes,collapse="_")
            grp[j]=paste(inHashes,collapse="_")
        }
        cd45Clust2=grp

        ## Mark cluster splits
        cellId=1:nrow(clustInfoThis)
        #z1=cutree(clustObjThis,k=80)
        #z1=rep("0",length(cd45Clust)); z1[cellId]=subClust
        cellId=which(cd45Clust%in%snacsObj$annHash$hashNames)
        z1=rep("0",length(cd45Clust)); z1[cellId]=x
        j=clustObjThis$order
        grp=as.integer(as.factor(z1[j]))
        j=which(diff(grp)!=0)
        jj=j
        j1=c(1,j+1); j2=c(j,length(grp))
        for (j in 1:length(j1)) {
            if (any(grp[j1[j]:j2[j]]!=grp[j1[j]])) print(j)
        }
        clustSplitRnd2=rep(1,length(grp))
        for (j in 1:length(j1)) {
            if (j%%2==0) clustSplitRnd2[j1[j]:j2[j]]=2
        }
        clustSplitRnd2=clustSplitRnd2[match(clustInfo$id,clustObjThis$label[clustObjThis$order])]
        subClustRnd2=z1
    }

    ###########################################################
    ###########################################################
    
    snacsObj$annCell$hashCall=""
    j=match(snacsObj$hclustObj_bestSNPs$label,snacsObj$annCell$id); j1=which(!is.na(j)); j2=j[j1]
    snacsObj$annCell$hashCall[j2]=cd45Clust[j1]
    snacsObj$annCell$hashCall[snacsObj$hashCall==""]=NA
    snacsObj$annCell$clustSplit=clustSplit
    snacsObj$annCell$subClust=subClust
    snacsObj$centroid=centroid
    snacsObj$dist2centroidMat=dist2centroidMat
    if (is.finite(minClustSizeRnd2)) {
        snacsObj$annCell$hashCallRnd2[j2]=cd45Clust2[j1]
        snacsObj$annCell$hashCallRnd2[snacsObj$hashCallRnd2==""]=NA
        snacsObj$annCell$clustSplitRnd2=clustSplitRnd2
        snacsObj$annCell$subClustRnd2=subClustRnd2
    }
    #print(table(duplicated(snacsObj$annCell$subClustRnd2)))
    
    cat('"hashCall" column added to "annCell" table in SNACS object\n',sep="")
    
    cat("\n")
    print(table(cd45Clust,cluster=clustInfo[match(snacsObj$hclustObj_bestSNPs$label,clustInfo$id),"clustId"],exclude=NULL,dnn=list("Hash calls",paste0("Cell clusters with best SNPs"))))

    invisible(snacsObj)
}

####################################################################
####################################################################
#' Generate hash antibody background density plot.
#'
#' The density plots are helpful is determining the background distribution of a hash antibody. The plots are saved in "../output" folder of a "pdf" or "png" format is specified. There are 3 figures generated for each subject. The top figure shows the distribution of a hash antibody. The middle figure shows the distribution of the background only. The red lines in the top and middle figures show the distribution of the estimated background of the hash antibody. The green lines mark the median and 95th percentile of the background distribution. The bottom figure is the histogram of the hash antibody data. The distribution of the antibody measure of a hash is expected to be a mixture of a background and a foreground signal where the latter signals the presence of that hash antibody. The background of the hash antibody is estimated by determining the background mode based on all cells and generating an empirical symmetric distribution around it.

#'
#' @param snacsObj SNACSList object
#' @param backgndThreshold Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell. Default is 0.95. Range is 0-1
#' @param cellProportionBelowBackgndMode Maximum proportion of cells which can be below the mode of the estimated hash background distribution. Default is 0.6. Range is 0-1
#' @param cellProportionForModeDetection Proportion of cells used to estimate mode of the background distribution. Used only if "cellProportionBelowBackgndMode" threshold is not met; otherwise, all cells are used. Default is 0.75. Range is 0-1
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @export
generateHashDensityPlot=function(snacsObj,backgndThreshold=0.95,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,outputFormat=c("","pdf","png")) {
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
    clustInfo=data.frame(id=snacsObj$annCell$id,t(snacsObj$hashes),stringsAsFactors=F)
    clustInfo=clustInfo[match(snacsObj$hclustObj_bestSNPs$labels,clustInfo$id),]

    ## --------------------------------------
    ## For manual annotation of doublets
    nClustFilt=NA; clustFilt=NA; clustRename=NA
    ## --------------------------------------

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
        if (stats::quantile(hashMat[,k],probs=cellProportionBelowBackgndMode)<x$x[k1]) k1=which.max(x$y[1:stats::quantile(1:(k1-1),probs=cellProportionForModeDetection)])
        xx=x$x[k1]
        x2=x$x[which(x$x<xx)]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBGSymData=x2
        hashBackgndECDF=stats::ecdf(hashBGSymData)
        x2=hashMat[which(hashMat[,k]<xx),k]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBGData=x2
        ecdfbg=1-hashBackgndECDF(hashBGData)
        vertLine=c(stats::median(x2),stats::quantile(x2,probs=backgndThreshold)); vertLineLab=c("med",backgndThreshold)
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
