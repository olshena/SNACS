###########################################################
###########################################################
#' Select SNPs that best separate the hashes
#'
#' Group cells into constituent samples based on hash antibody data. Select the best SNPs that separate the cell groups using the mutation data.
#'
#' @param snacsObj SNACSList object
#' @param numSNP Integer. Number of pairs of SNPs (high in one sample / low in second sample and vice versa) that best separate sample-pairs. There will be s times n(n-1) number of SNPs where s is the number of SNPs and n is the number of samples. Default is 3
#' @param backgndThreshold Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell. Default is 0.95. Range is 0-1
#' @param cellProportionAboveBackgnd Numeric. Proportion of cells in a cluster that have to have signal above the background of a hash, determined by "backgndThreshold" parameter, for the cluster to be assigned to the hash. Default is 0.5. Range is 0-1
#' @param cellProportionBelowBackgndMode Numeric. Maximum proportion of cells which can be below the mode of the estimated hash background distribution. Default is 0.6. Range is 0-1
#' @param cellProportionForModeDetection Numeric. Proportion of cells used to estimate mode of the background distribution. Used only if "cellProportionBelowBackgndMode" threshold is not met; otherwise, all cells are used. Default is 0.75. Range is 0-1
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @return A SNACSList object
#' @export
getBestSNPs=function(snacsObj,numSNP=3,backgndThreshold=0.95,cellProportionAboveBackgnd=0.5,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,outputFormat=c("","pdf","png")) {
    ## -----------------------------------
    if (is.na(match("filtered",snacsObj[["processLevels"]]))) stop("Run filterData() before selecting best SNPs\n")

    ## -----------------------------------
    snacsObj=clusterSampleWithAntibodyData(snacsObj,backgndThreshold=backgndThreshold,cellProportionBelowBackgndMode=cellProportionBelowBackgndMode,cellProportionForModeDetection=cellProportionForModeDetection)
    generateAntibodyDensityPlot(snacsObj,backgndThreshold=backgndThreshold,cellProportionBelowBackgndMode=cellProportionBelowBackgndMode,cellProportionForModeDetection=cellProportionForModeDetection,outputFormat=outputFormat)

    ## Association of mutation status of each SNP with each sample cluster pair
    samUniq=snacsObj$annHash$hashNames
    n=nrow(snacsObj$annSNP)*length(samUniq)*(length(samUniq)-1)/2
    tmp=rep(NA,n); tmpC=rep("",n)
    stat=data.frame(id=tmpC,sam1=tmpC,sam2=tmpC,stat=tmp,pv=tmp,propMutIn1stSample=tmp)
    k=1
    for (samId1 in 1:(length(samUniq)-1)) {
        for (samId2 in (samId1+1):length(samUniq)) {
            j=which(snacsObj$annCell$hashClust%in%c(samUniq[samId1],samUniq[samId2]))
            for (sId in 1:nrow(snacsObj$annSNP)) {
                stat$id[k]=snacsObj$annSNP$id[sId]
                stat$sam1[k]=samUniq[samId1]
                stat$sam2[k]=samUniq[samId2]
                x=table(snacsObj$annCell$hashClust[j],snacsObj$mut[sId,j])
                if (nrow(x)==2 & ncol(x)==2) {
                    fit=try(stats::chisq.test(x))
                    if (inherits(fit,"htest")) {
                        stat$stat[k]=fit$statistic
                        stat$pv[k]=fit$p.value
                        stat$propMutIn1stSample[k]=fit$observed[1,2]/sum(fit$observed[1,])
                    }
                }
                k=k+1
            }
        }
    }
    stat=stat[order(stat$pv,decreasing=F),]

    ## Top SNPs associated with sample cluster pairs
    n=2*numSNP*length(samUniq)*(length(samUniq)-1)/2
    tmp=rep(NA,n); tmpC=rep("",n)
    statBest=data.frame(id=tmpC,sam1=tmpC,sam2=tmpC,dirn=tmpC)
    k=0
    for (samId1 in 1:(length(samUniq)-1)) {
        for (samId2 in (samId1+1):length(samUniq)) {
            i=which(stat$sam1==samUniq[samId1] & stat$sam2==samUniq[samId2] & stat$propMutIn1stSample>0.5)[1:numSNP]
            k2=k+(1:numSNP)
            statBest$id[k2]=stat$id[i]
            statBest$sam1[k2]=samUniq[samId1]
            statBest$sam2[k2]=samUniq[samId2]
            statBest$dirn[k2]=paste0(samUniq[samId1],"up_",samUniq[samId2],"down")
            i=which(stat$sam1==samUniq[samId1] & stat$sam2==samUniq[samId2] & stat$propMutIn1stSample<0.5)[1:numSNP]
            k=k+numSNP
            k2=k+(1:numSNP)
            statBest$id[k2]=stat$id[i]
            statBest$sam1[k2]=samUniq[samId1]
            statBest$sam2[k2]=samUniq[samId2]
            statBest$dirn[k2]=paste0(samUniq[samId2],"up_",samUniq[samId1],"down")
            k=k+numSNP
        }
    }
    table(uniqueSNP=!duplicated(statBest$id))
    id=unique(statBest$id[duplicated(statBest$id)])
    if (length(id)!=0) {
        for (ii in 1:length(id)) {
            i=which(statBest$id==id[ii])
            statBest$dirn[i]=unique(paste0(statBest$dirn[i],collapse="|"))
        }
    }
    statBest=statBest[!duplicated(statBest$id),]

    i=match(statBest$id,snacsObj$annSNP$id)
    snacsObj$mut=snacsObj$mut[i,]; snacsObj$annSNP=snacsObj$annSNP[i,]
    snacsObj$annSNP$dirn=statBest$dirn
    snacsObj[["processLevels"]]=c(snacsObj[["processLevels"]],"bestSNPs")
    if (!is.null(snacsObj$depthTotal)) {
        snacsObj$depthTotal=snacsObj$depthTotal[i,]
        snacsObj$depthAlt=snacsObj$depthAlt[i,]
    }

    invisible(snacsObj)
}

###########################################################
###########################################################
#' Cluster cells based on the SNP data
#'
#' Cluster cells, where the number of clusters is the number of constituent samples, using the mutation data with the selected SNPs. Generate heatmap of the mutation data.
#'
#' @param snacsObj SNACSList object
#' @param clustMethod Character. Clustering method for ranking SNPs. Options are "hclust" and "skmean". Default is "hclust"
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @return A SNACSList object
#' @export
clusterCellsWithSNPdata=function(snacsObj,clustMethod=c("hclust","skmean"),outputFormat=c("","pdf","png")) {
    clustMethod=clustMethod[1]
    outputFormat=outputFormat[1]

    ## -----------------------------------
    if (is.na(match("bestSNPs",snacsObj[["processLevels"]]))) stop("Run getBestSNPs() before clustering SNPs\n")
    if (is.na(match("imputed",snacsObj[["processLevels"]]))) stop("Run imputeMissingMutations() before clustering SNPs\n")

    pvSnpThres=NA
    cellClusterFileName=NA
    subsetSnpFlag=""
    subsetCellFlag=""
    
    outputFileName=paste0("heatmap",ifelse(clustMethod=="skmean","_skmean",""),"_pvRanked_",snacsObj$exptName)
    if (outputFormat=="") h_title="Heatmap of ranked SNPs" else h_title=NULL
    
    cell_anno_var_rankedSNPs=snacsObj$annHash$hashNames
    
    ################################################
    ## Compute cosine cross-distances between the rows of matrices
    getSKmeansDist=function(x) {stats::as.dist(skmeans::skmeans_xdist(x))}
    distfun=getSKmeansDist
    linkMethod="ward.D2"

    ################################################
    datThis=snacsObj$mut
    annSNPthis=snacsObj$annSNP
    annCellThis=snacsObj$annCell
    hashesThis=snacsObj$hashes

    ################################################
    if (is.na(pvSnpThres)) {
        if (clustMethod=="skmean") {
            clustCell=skmeans::skmeans(t(datThis),k=nrow(snacsObj$annHash))
            silCell=cluster::silhouette(clustCell)
            clustCell=clustCell$cluster[match(annCellThis$id,names(clustCell$cluster))]
            silCell=silCell[match(annCellThis$id,rownames(silCell)),]
            annCellThis$cluster_skmean=clustCell
            annCellThis$cluster_silhouette=silCell[,"sil_width"]
            grp=annCellThis$cluster_skmean
        } else if (!is.na(cellClusterFileName)) {
            clustCell=utils::read.table(cellClusterFileName, sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=-1)
            j=match(annCellThis$id,clustCell$id)
            if (any(is.na(j))) stop("filter_data_for_heatmap: Cell mismatch !!!")
            clustCell=clustCell[j,]
            annCellThis$cluster_hclust=clustCell[,paste0("clustId_",nrow(snacsObj$annHash))]
            grp=annCellThis$cluster_hclust
        } else {
            distMat=distfun(t(datThis))
            clustCell=stats::hclust(distMat,method=linkMethod)
            clustCell=heatmap4::cutCluster(clustCell,ann=annCellThis,nClust=nrow(snacsObj$annHash))
            j=match(annCellThis$id,clustCell$id)
            if (any(is.na(j))) stop("filter_data_for_heatmap: Cell mismatch !!!")
            clustCell=clustCell[j,]
            annCellThis$cluster_hclust=clustCell[,paste0("clustId_",nrow(snacsObj$annHash))]
            grp=annCellThis$cluster_hclust
        }
        grpUniq=unique(grp)
        pvMat=matrix(nrow=nrow(annSNPthis),ncol=length(grpUniq))
        lorMat=matrix(nrow=nrow(annSNPthis),ncol=length(grpUniq))
        for (gId in 1:length(grpUniq)) {
            grp2=rep(1,length(grp))
            grp2[which(!grp%in%grpUniq[gId])]=0
            for (i in 1:nrow(annSNPthis)) {
                x=table(datThis[i,],grp2)
                if (nrow(x)<2 | ncol(x)<2) {
                    pvMat[i,gId]=99
                } else {
                    res=stats::fisher.test(datThis[i,],grp2)
                    pvMat[i,gId]=res$p.value
                    lorMat[i,gId]=log(res$estimate)
                }
            }
        }
        annSNPthis$minPValue_1VsOtherClusters=apply(pvMat,1,min,na.rm=T)
        annSNPthis$maxLOR_1VsOtherClusters=apply(lorMat,1,max,na.rm=T)
        i=rev(order(annSNPthis$minPValue_1VsOtherClusters))
    } else {
        i=which(annSNPthis$minPValue_1VsOtherClusters<=pvSnpThres)
        if (length(i)==0) {
            cat("\n\n No. of SNPs :",nrow(annSNPthis),"\n",sep="")
            cat("SNP p-values:\n")
            print(summary(annSNPthis$minPValue_1VsOtherClusters))
            stop(paste0("No SNPs with p-value <= ",pvSnpThres," !!!"))
        }
    }
    datThis=datThis[i,]
    annSNPthis=annSNPthis[i,]
    if (!is.null(snacsObj$depthTotal)) {
        snacsObj$depthTotal=snacsObj$depthTotal[i,]
        snacsObj$depthAlt=snacsObj$depthAlt[i,]
    }

    j=apply(datThis,2,function(x) {mean(x==1)}); j=which(j>0 & j<1) ## Exclude cells with no or all mutations
    datThis=datThis[,j]
    hashesThis=hashesThis[,j]
    annCellThis=annCellThis[j,]
    if (!is.null(snacsObj$depthTotal)) {
        snacsObj$depthTotal=snacsObj$depthTotal[,j]
        snacsObj$depthAlt=snacsObj$depthAlt[,j]
    }

    if (clustMethod=="skmean") {
        if (!is.na(pvSnpThres)) {
            clustCell=skmeans::skmeans(t(datThis),k=nrow(snacsObj$annHash))
            silCell=cluster::silhouette(clustCell)
            clustCell=clustCell$cluster[match(annCellThis$id,names(clustCell$cluster))]
            silCell=silCell[match(annCellThis$id,rownames(silCell)),]
            annCellThis$cluster_skmean2=clustCell
            annCellThis$cluster_silhouette2=silCell[,"sil_width"]
            grp=annCellThis$cluster_skmean2
        }
        grp=annCellThis$cluster_skmean
        grpUniq=unique(grp)
        j=c()
        for (gId in 1:length(grpUniq)) {
            jj=which(grp==grpUniq[gId])
            k=order(annCellThis$cluster_silhouette[jj])
            j=c(j,jj[k])
        }
        datThis=datThis[,j]
        hashesThis=hashesThis[,j]
        annCellThis=annCellThis[j,]
        if (!is.null(snacsObj$depthTotal)) {
            snacsObj$depthTotal=snacsObj$depthTotal[,j]
            snacsObj$depthAlt=snacsObj$depthAlt[,j]
        }
    }

    ################################################
    if (subsetSnpFlag!="") {
        if (subsetSnpFlag=="_knownHashSnpOrder") {
            tbl=utils::read.table(paste0("../heatmap/_knownHash/clustInfoSnp_knownHash_",snacsObj$exptName,".txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            i=match(tbl$id,annSNPthis$id)
            i=rev(i)
            datThis=datThis[i,]
            annSNPthis=annSNPthis[i,]
            if (!is.null(snacsObj$depthTotal)) {
                snacsObj$depthTotal=snacsObj$depthTotal[i,]
                snacsObj$depthAlt=snacsObj$depthAlt[i,]
            }
        }
    }
    
    ################################################
    j=apply(datThis,2,function(x) {mean(x==1)}); j=which(j>0 & j<1) ## Exclude cells with no or all mutations
    datThis=datThis[,j]
    hashesThis=hashesThis[,j]
    annCellThis=annCellThis[j,]
    if (!is.null(snacsObj$depthTotal)) {
        snacsObj$depthTotal=snacsObj$depthTotal[,j]
        snacsObj$depthAlt=snacsObj$depthAlt[,j]
    }

    ################################################
    names(annCellThis)[match(paste0("cluster_",clustMethod),names(annCellThis))]=paste0("clustRankedSNPs_",clustMethod)

    snacsObj[["mut"]]=datThis
    snacsObj[["hashes"]]=hashesThis
    snacsObj[["annCell"]]=annCellThis
    snacsObj[["annSNP"]]=annSNPthis

    ################################################
    clustObj=createHeatmap(snacsObj,cell_anno_var=cell_anno_var_rankedSNPs,col_dend=T,row_dend=F,h_title=h_title,outputFormat=outputFormat,outputFileName=outputFileName)
    snacsObj[["hclustObj_bestSNPs"]]=clustObj$colClust

    ################################################
    invisible(snacsObj)
}

####################################################################
####################################################################
#' Make hash calls
#'
#' Make hash calls based on mutation and hash antibody data. The cluster information from the mutation data and the hash antibody data are used to make hash calls. Each cell is assigned one of the hash IDs (a single sample) or a set of hash IDs (a multiplet), when it is a mixture of samples. The distribution of the antibody measure of a hash is expected to be a mixture of a background and a foreground signal where the latter signals the presence of antibody for the hash. The clusters based on mutation data give a rough assignment of hash calls. The hash calls are refined by looking at the antibody data of the cells which cluster together. The cell clusters obtained from running "clusterCellsWithSNPdata" are split into sub-clusters and each sub-cluster is assigned a hash ID if enough cells in that group have a signal above the background for the corresponding hash. The sub-clusters are created using the cell clusters from the second round of hierarchical clustering of the mutation data. First the antibody measures of the top two cluster pairs are compared using Hotelling's T2 test. If there is significant difference , the two clusters are splits into sub-clusters and each of those cluster pairs are tested for difference. The tree is traversed until there are no pairs of significantly different cluster pairs. The splitting of a cluster stops when it reaches a minimum number of cells. The singlet and multiplet assignments are made to each of the sub-clusters using the hash antibody data. The antibody expression for a sample is expected to have a bimodal distribution representing a background, which are cells not from the sample, and a foreground, which are cells belonging to the sample. The background distribution is estimated by considering the data from the lower bound till the mode of the background and generating a symmetrical distribution from it. Then, the empirical cumulative distribution function of this distribution becomes the estimated background distribution. A sub-cluster is assigned the ID of each sample which has at least a certain proportion of cells, defined by the parameter "cellProportionAboveBackgnd", with antibody expression above a certain quantile, defined by the parameter "backgndThreshold", of the hash background. A sub-cluster with multiple samples above their corresponding hash background is considered a multiplet. The result is stored in the column "hashCallRnd1" in the "annCell" table of the SNACS object.
#' The above hash calls are refined by detecting narrow regions of multiplets. The cells, ordered as in the clusters above, are segmented by applying circular binary segmentation (CBS) on the distance to the centroid of the hash call groups. A segment is assigned the ID of each hash which has at least a certain proportion of cells, defined by the parameter "cellProportionAboveBackgnd", with antibody expression above a certain quantile, defined by the parameter "backgndThresRnd2", of the hash background. A segment with multiple hashes above their corresponding hash background is considered a multiplet. Only narrow segments, specified by "maxClustSizeRnd2", are considered. The result is stored in the column "hashCallRnd2" in the "annCell" table of the SNACS object
#'
#' @param snacsObj SNACSList object
#' @param backgndThreshold Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell. Default is 0.95. Range is 0-1
#' @param cellProportionAboveBackgnd Numeric. Proportion of cells in a cluster that have to have signal above the background of a hash, determined by "backgndThreshold" parameter, for the cluster to be assigned to the hash. Default is 0.5. Range is 0-1
#' @param cellProportionBelowBackgndMode Numeric. Maximum proportion of cells which can be below the mode of the estimated hash background distribution. Default is 0.6. Range is 0-1
#' @param cellProportionForModeDetection Numeric. Proportion of cells used to estimate mode of the background distribution. Used only if "cellProportionBelowBackgndMode" threshold is not met; otherwise, all cells are used. Default is 0.75. Range is 0-1
#' @param minClustSize Numeric. Minimum number of cells required to be in a cluster. Default is 2
#' @param clustComparePValue Numeric. P-value threshold to compare cluster pairs. Default is 10^-5
#' @param maxClustSampleSize Numeric. Maximum number of cells in a cluster that will be used to compare. If more, then this number of cells will be sampled. Default is Inf
#' @param clustCompareMethod Character. Test used to compare clusters. Default is t-test
#' @param makeHashCallRnd2 Logical. Make second round of hash call to detect narrow multiplet regions. Default is TRUE
#' @param maxClustSizeRnd2 Integer. Maximum number of cells required to be in a cluster for making second round of hash calls. Default is 100
#' @param backgndThresRnd2 Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell for making second round of hash calls. Default is 0.75. Range is 0-1
#' @param dataTypeRnd2 Character. Type of data to be used for splitting the cells when making second round of hash calls. Default is euclidean distance of the cells to the cluster means from first round of hash calls
#' @param cbsAlpha Numeric. Significance level passed to DNAcopy::segment to for splitting the cells when making second round of hash calls. Default is 0.1
#' @return A SNACSList object
#' @export
makeHashCall=function(snacsObj,backgndThreshold=0.95,cellProportionAboveBackgnd=0.5,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,minClustSize=2,clustComparePValue=10^-5,maxClustSampleSize=Inf,clustCompareMethod=c("t","hotelling"),makeHashCallRnd2=TRUE,maxClustSizeRnd2=100,backgndThresRnd2=0.75,dataTypeRnd2=c("euclidean","sum of squares","log2 euclidean","log2 sum of squares"),cbsAlpha=0.1) {
    clustCompareMethod=clustCompareMethod[1]
    dataTypeRnd2=dataTypeRnd2[1]

    ## -----------------------------------
    if (is.na(match("hclustObj_bestSNPs",names(snacsObj)))) stop("Run clusterCellsWithSNPdata() before making hash calls\n")

    ## -----------------------------------
    clustInfo=data.frame(id=snacsObj$annCell$id,t(snacsObj$hashes),stringsAsFactors=F)
    clustInfo=clustInfo[match(snacsObj$hclustObj_bestSNPs$labels,clustInfo$id),]

    ## --------------------------------------
    clustInfo$clustId=stats::cutree(snacsObj$hclustObj_bestSNPs,k=nrow(snacsObj$annHash))
    
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
    
    ###########################################################
    ## Make hash call round 1
    
    clustInfo$hashCallRnd1="Multiplet"
    grp=clustInfo$clustId
    grpUniq=sort(unique(grp))
    for (gId in 1:length(grpUniq)) {
        j=which(grp==grpUniq[gId])
        x=apply(hashMat[j,],2,mean,na.rm=T)
        clustInfo$hashCallRnd1[j]=names(x)[which.max(x)]
    }

    hashCallRnd1=clustInfo$hashCallRnd1
    
    clustObjThis=snacsObj$hclustObj_bestSNPs
    j=match(clustObjThis$label,clustInfo$id)
    clustInfoThis=clustInfo[j,]
    hashMatThis=hashMat[j,]

    heightUniq=rev(unique(clustObjThis$height))
    cellId=1:nrow(clustInfoThis)
    value=rep(NA,length(cellId))
    valueOrig=rep(NA,length(cellId))
    grpHigher=stats::cutree(clustObjThis,k=1); grpHigher=grpHigher[cellId]
    hMaxId=max(2,length(heightUniq)-1)
    for (hId in 2:hMaxId) {
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
                        res=DescTools::HotellingsT2Test(hashMatThis[j,]~grpThis)
                        pv=res$p.value
                    } else {
                        pv=rep(NA,ncol(hashMatThis))
                        for (k in 1:ncol(hashMatThis)) {
                            res=stats::t.test(hashMatThis[j,k]~grpThis,var.equal=T)
                            pv[k]=res$p.value
                        }
                        pv=min(pv*ncol(hashMatThis))
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
        inSamples=c()
        for (sampleId in snacsObj$annHash$hashNames) {
            if (round(mean(hashMatThis[j,sampleId]>hashBackgnd["thres",sampleId]),2)>=cellProportionAboveBackgnd) {inSamples=c(inSamples,sampleId)}
        }
        hashCallRnd1[j]=paste(inSamples,collapse="_")
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
    clustSplitRnd1=rep(1,length(grp))
    for (j in 1:length(j1)) {
        if (j%%2==0) clustSplitRnd1[j1[j]:j2[j]]=2
    }
    clustSplitRnd1=clustSplitRnd1[match(clustInfo$id,clustObjThis$label[clustObjThis$order])]
    subClustRnd1=x

    ## Distance to centroid
    grp=hashCallRnd1
    grpUniq=sort(unique(grp))
    centroid=matrix(nrow=ncol(hashMatThis),ncol=length(grpUniq),dimnames=list(colnames(hashMatThis),grpUniq))
    for (gId in 1:length(grpUniq)) {
        j=which(grp==grpUniq[gId])
        centroid[,gId]=apply(hashMatThis[j,],2,mean,na.rm=T)
    }
    dist2centroidMat=matrix(nrow=length(grpUniq),ncol=nrow(hashMatThis),dimnames=list(grpUniq,rownames(hashMatThis)))
    j=which(grp%in%colnames(hashMatThis))
    j=1:nrow(hashMatThis)
    for (gId in 1:length(grpUniq)) {
        dist2centroidMat[gId,j]=apply(hashMatThis[j,],1,function(x,centroid) {sum((x-centroid)^2)},centroid=centroid[,gId])
    }
    switch(dataTypeRnd2,
        "euclidean"={dist2centroidMat=sqrt(dist2centroidMat)},
        "log2 sum of squares"={dist2centroidMat=log2(dist2centroidMat)},
        "log2 euclidean"={dist2centroidMat=log2(sqrt(dist2centroidMat))}
    )

    ###########################################################
    ## Make hash call round 2
    
    if (makeHashCallRnd2) {
        ## Detect and assign narrow regions as multiplets
        ## Make splits based on cbs on distance to centroids
        
        cellId=which(hashCallRnd1%in%snacsObj$annHash$hashNames)
        #x1=hashMatThis[clustObjThis$order,]
        x1=t(dist2centroidMat[,clustObjThis$order])
        x2=DNAcopy::CNA(genomdat=x1,
        chrom=rep(1,nrow(x1)),maploc=1:nrow(x1),
                          data.type="logratio",sampleid=colnames(x1))
        x3=DNAcopy::segment(x2,alpha=cbsAlpha,verbose=0)
        x3=x3$output
        x=rep(NA,nrow(x1))
        names(x)=clustObjThis$label[clustObjThis$order]
        grp=hashCallRnd1[clustObjThis$order]
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

        grp=hashCallRnd1
        #x[is.na(x)]="0"
        for (k in unique(x[!is.na(x)])) {
            j=cellId[which(x==k)]
            if (length(j)>maxClustSizeRnd2) next
            inSamples=c()
            for (sampleId in snacsObj$annHash$hashNames) {
                if (round(mean(hashMatThis[j,sampleId]>hashBackgnd["thresRnd2",sampleId]),2)>=cellProportionAboveBackgnd) {inSamples=c(inSamples,sampleId)}
            }
            grp[j]=paste(inSamples,collapse="_")
        }
        hashCallRnd2=grp

        ## Mark cluster splits
        cellId=1:nrow(clustInfoThis)
        cellId=which(hashCallRnd1%in%snacsObj$annHash$hashNames)
        z1=rep("0",length(hashCallRnd1)); z1[cellId]=x
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
    } ## end of makeHashCallRnd2

    ###########################################################
    ###########################################################
    
    snacsObj$annCell$hashCallRnd1=""
    j=match(snacsObj$hclustObj_bestSNPs$label,snacsObj$annCell$id); j1=which(!is.na(j)); j2=j[j1]
    snacsObj$annCell$hashCallRnd1[j2]=hashCallRnd1[j1]
    snacsObj$annCell$hashCallRnd1[snacsObj$hashCallRnd1==""]=NA
    snacsObj$annCell$clustSplitRnd1=clustSplitRnd1
    snacsObj$annCell$subClustRnd1=subClustRnd1
    snacsObj$centroid=centroid
    snacsObj$dist2centroidMat=dist2centroidMat
    if (makeHashCallRnd2) {
        snacsObj$annCell$hashCallRnd2=""
        snacsObj$annCell$hashCallRnd2[j2]=hashCallRnd2[j1]
        snacsObj$annCell$hashCallRnd2[snacsObj$hashCallRnd2==""]=NA
        snacsObj$annCell$clustSplitRnd2=clustSplitRnd2
        snacsObj$annCell$subClustRnd2=subClustRnd2
    }
    
    cat('"hashCall" column(s) added to "annCell" table in SNACS object\n',sep="")
    
    cat("\n")
    if (makeHashCallRnd2) {hashCall=hashCallRnd2
    } else {hashCall=hashCallRnd1}
    print(table(hashCall,cluster=clustInfo[match(snacsObj$hclustObj_bestSNPs$label,clustInfo$id),"clustId"],exclude=NULL,dnn=list("hash calls",paste0("Cell clusters with best SNPs"))))
    rm(hashCall)

    invisible(snacsObj)
}

####################################################################
####################################################################
#' Assign sample IDs to cells based on hash antibody data
#'
#' Assign sample IDs to cells based on hash antibody data
#'
#' @param snacsObj SNACSList object
#' @param backgndThreshold Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell. Default is 0.95. Range is 0-1
#' @param cellProportionBelowBackgndMode Numeric. Maximum proportion of cells which can be below the mode of the estimated hash background distribution. Default is 0.6. Range is 0-1
#' @param cellProportionForModeDetection Numeric. Proportion of cells used to estimate mode of the background distribution. Used only if "cellProportionBelowBackgndMode" threshold is not met; otherwise, all cells are used. Default is 0.75. Range is 0-1
#' @param hashThreshold Numeric. Threshold of the hash value above which the antibody will be considered to be expressed in a cell. Default is 0.5. Used if backgndThreshold = NA
#' @export
clusterSampleWithAntibodyData=function(snacsObj,backgndThreshold=0.95,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,hashThreshold=0.5) {
    ## --------------------------------------
    dirOutput="../output/hashPairPlot/"
    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))

    ## --------------------------------------
    ## --------------------------------------
    backgndThresRnd2=0.75
    
    hashMat=t(snacsObj$hashes)
    ## --------------------------------------
    ## --------------------------------------

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
    
    hashCallThis=rep("",nrow(snacsObj$annCell))
    hashMatThis=hashMat
    j=1:length(hashCallThis)
    for (sampleId in snacsObj$annHash$hashNames) {
        if (is.na(backgndThreshold)) {
            jj=which(hashMatThis[,sampleId]>=hashThreshold)
        } else {
            jj=which(hashMatThis[,sampleId]>hashBackgnd["thres",sampleId])
        }
        if (length(jj)!=0) hashCallThis[jj]=paste0(hashCallThis[jj],"_",sampleId)
    }
    hashCallThis=sub("^_", "",hashCallThis)
    
    snacsObj$annCell$hashClust=hashCallThis

    invisible(snacsObj)
}

####################################################################
####################################################################
#' Generate hash antibody background density plot
#'
#' The density plots are helpful is determining the background distribution of a hash antibody. The plots are saved in "../output" folder of a "pdf" or "png" format is specified. There are 3 figures generated for each subject. The top figure shows the distribution of a hash antibody. The middle figure shows the distribution of the background only. The red lines in the top and middle figures show the distribution of the estimated background of the hash antibody. The green lines mark the median and 95th percentile of the background distribution. The bottom figure is the histogram of the hash antibody data. The distribution of the antibody measure of a hash is expected to be a mixture of a background and a foreground signal where the latter signals the presence of that hash antibody. The background of the hash antibody is estimated by determining the background mode based on all cells and generating an empirical symmetric distribution around it.

#'
#' @param snacsObj SNACSList object
#' @param backgndThreshold Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell. Default is 0.95. Range is 0-1
#' @param cellProportionBelowBackgndMode Numeric. Maximum proportion of cells which can be below the mode of the estimated hash background distribution. Default is 0.6. Range is 0-1
#' @param cellProportionForModeDetection Numeric. Proportion of cells used to estimate mode of the background distribution. Used only if "cellProportionBelowBackgndMode" threshold is not met; otherwise, all cells are used. Default is 0.75. Range is 0-1
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @export
generateAntibodyDensityPlot=function(snacsObj,backgndThreshold=0.95,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,outputFormat=c("","pdf","png")) {
    outputFormat=outputFormat[1]

    ## -----------------------------------
    if (outputFormat=="") {
        plotInfo=list(cexMain=1,cexLab=1,cexAxis=1)
    } else {
        dirResult="../output/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
        dirResult="../output/hashDensityPlot/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
        plotInfo=list(cexMain=2,cexLab=1.5,cexAxis=1.5)
    }
    
    ## --------------------------------------
    dirOutput="../output/hashPairPlot/"
    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))
    
    hashMat=t(snacsObj$hashes)
    ## --------------------------------------
    ## --------------------------------------

    hashBackgnd=matrix(nrow=2,ncol=nrow(snacsObj$annHash),dimnames=list(c("mean","sd"),snacsObj$annHash$hashNames))
    ylimHist=rep(NA,ncol(hashMat))
    for (k in 1:ncol(hashMat)) {
        y=graphics::hist(hashMat[,k],breaks=diff(range(hashMat[,k]))/.1,plot=F)
        ylimHist[k]=max(y$counts)
    }
    ylimHist=c(0,max(ylimHist))
    for (k in 1:ncol(hashMat)) {
        switch(outputFormat,
            "png"={grDevices::png(paste0(dirResult,"densityPlot_hash_",snacsObj$annHash$hashNames[k],"_",snacsObj$exptName,".png"),width=5*240,height=3*240)},
            "pdf"={grDevices::pdf(paste0(dirResult,"densityPlot_hash_",snacsObj$annHash$hashNames[k],"_",snacsObj$exptName,".pdf"),width=7,height=7)}
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
        graphics::plot(x,xlim=xlim,ylim=ylim,main=paste0(snacsObj$exptName,": ",colnames(hashMat)[k]),xlab="Hash",cex.main=plotInfo$cexMain,cex.lab=plotInfo$cexLab,cex.axis=plotInfo$cexAxis)
        graphics::lines(stats::density(x2,bw="SJ"),col="red"); graphics::abline(v=vertLine,col="green")
        graphics::plot(stats::density(x2,bw="SJ"),xlim=xlim,ylim=ylim,main="",xlab="Hash background",col="red",cex.main=plotInfo$cexMain,cex.lab=plotInfo$cexLab,cex.axis=plotInfo$cexAxis); graphics::abline(v=vertLine,col="green")
        graphics::axis(side=3,at=vertLine,labels=vertLineLab,cex.axis=plotInfo$cexAxis,las=3,col="green")
        hashBackgnd[,k]=c(mean(x2),stats::sd(x2))
        x=hashMat[,k]
        graphics::hist(x,freq=T,breaks=diff(range(x))/.1,xlim=xlim,ylim=ylimHist,main=paste0(snacsObj$exptName,": ",colnames(hashMat)[k]),xlab="Hash",ylab="Count",cex.main=plotInfo$cexMain,cex.lab=plotInfo$cexLab,cex.axis=plotInfo$cexAxis)

        if (outputFormat!="") grDevices::dev.off()
    }
}

###########################################################
###########################################################
