#' Run SNACS to make demultiplex samples in a multiplexed SNP experiment
#'
#' Make calls for constituent samples in a multiplexed SNP experiment using SNP and hash antibody data.
#'
#' @param snacsObj SNACSList object
#' @param proportionMissingPerSNP Numeric. Only SNPs with lower than this proportion of missing values are considered. Default is 0.4. Range is 0-1
#' @param proportionMissingPerCell Numeric. Only cells with lower than this proportion of missing values are considered. Default is 0.4. Range is 0-1
#' @param proportionMutatedPerSNP Numeric vector. Only SNPs having proportion of mutations within this range are considered. Default is c(0.05,0.95). Range is 0-1
#' @param proportionMutatedPerCell Numeric vector. Only cells having proportion of mutations within this range are considered. Default is c(0,1). Range is 0-1
#' @param numSNP Integer. Number of pairs of SNPs (high in one sample / low in second sample and vice versa) that best separate sample-pairs. There will be s times n(n-1) number of SNPs where s is the number of SNPs and n is the number of samples. Default is 3
#' @param minSampleSize Integer. Minimum number of cells that is estimated by SNACS to be in a constituent sample. Gives a warning if below this number. Default is 100
#' @param bgndThresDetMethod Character. Method for detecting threshold of background antibody distribution. Default is "automatic"
#' @param backgndThreshold Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell for making first round of SNACS calls. Default is 0.95. Range is 0-1
#' @param backgndThresRnd2 Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell for making second round of SNACS calls. Default is 0.75. Range is 0-1
#' @param cellProportionAboveBackgnd Numeric. Proportion of cells in a cluster that have to have signal above the background of a hash, determined by "backgndThreshold" parameter, for the cluster to be assigned to the hash. Default is 0.5. Range is 0-1
#' @param cellProportionBelowBackgndMode Numeric. Maximum proportion of cells which can be below the mode of the estimated hash background distribution. Default is 0.6. Range is 0-1
#' @param cellProportionForModeDetection Numeric. Proportion of cells used to estimate mode of the background distribution. Used only if "cellProportionBelowBackgndMode" threshold is not met; otherwise, all cells are used. Default is 0.75. Range is 0-1
#' @param hashThreshold Numeric. Threshold of the hash value above which the antibody will be considered to be expressed in a cell. If NA then backgndThreshold is used. Default is 0.5
#' @param bgndQuantileThreshold Numeric. Threshold of the hash value above which the antibody will be considered to be expressed in a cell in first round of SNACS call. Default is NA. Used if bgndThresDetMethod = "manual"
#' @param bgndQuantileThresRnd2 Numeric. Threshold of the hash value above which the antibody will be considered to be expressed in a cell in second round of SNACS call. Default is NA. Used if bgndThresDetMethod = "manual"
#' @param minClustSize Numeric. Minimum number of cells required to be in a cluster. Default is 2
#' @param clustComparePValue Numeric. P-value threshold to compare cluster pairs. Default is 10^-5
#' @param maxClustSampleSize Numeric. Maximum number of cells in a cluster that will be used to compare. If more, then this number of cells will be sampled. Default is Inf
#' @param clustCompareMethod Character. Test used to compare clusters. Default is t-test
#' @param maxClustSizeRnd2 Integer. Maximum number of cells required to be in a cluster for making second round of SNACS calls. Default is 100
#' @param dataTypeRnd2 Character. Type of data to be used for splitting the cells when making second round of SNACS calls. Default is euclidean distance of the cells to the cluster means from first round of SNACS calls
#' @param cbsAlpha Numeric. Significance level passed to DNAcopy::segment to for splitting the cells when making second round of SNACS calls. Default is 0.1
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @param verbose Logical. Prints information when running the method. Default is FALSE
#' @return A SNACSList object
#' @export
runSNACS=function(snacsObj,proportionMissingPerSNP=0.4,proportionMissingPerCell=0.4,proportionMutatedPerSNP=c(0.05,0.95),proportionMutatedPerCell=c(0,1),numSNP=3,minSampleSize=100,bgndThresDetMethod=c("automatic","manual","default modes","two modes"),backgndThreshold=0.95,backgndThresRnd2=0.75,cellProportionAboveBackgnd=0.5,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,hashThreshold=0.5,bgndQuantileThreshold=NA,bgndQuantileThresRnd2=NA,minClustSize=2,clustComparePValue=10^-5,maxClustSampleSize=Inf,clustCompareMethod=c("t","hotelling"),maxClustSizeRnd2=100,dataTypeRnd2=c("euclidean","sum of squares","log2 euclidean","log2 sum of squares"),cbsAlpha=0.1,outputFormat=c("","pdf","png"),verbose=FALSE) {
    bgndThresDetMethod=bgndThresDetMethod[1]
    clustCompareMethod=clustCompareMethod[1]
    dataTypeRnd2=dataTypeRnd2[1]
    outputFormat=outputFormat[1]

    snacsObj=filterData(snacsObj,
        proportionMissingPerSNP=proportionMissingPerSNP,
        proportionMissingPerCell=proportionMissingPerCell,
        proportionMutatedPerSNP=proportionMutatedPerSNP,
        proportionMutatedPerCell=proportionMutatedPerCell,
        verbose=verbose)

    done=T
    bgndThresDetMethodVec=bgndThresDetMethod
    if (bgndThresDetMethod=="automatic") bgndThresDetMethodVec=c("automatic","two modes")
    for (bgndThresDetMethodThis in bgndThresDetMethodVec) {
        snacsObj=getBestSNPs(snacsObj,
            numSNP=numSNP,
            minSampleSize=minSampleSize,
            bgndThresDetMethod=bgndThresDetMethod,
            backgndThreshold=backgndThreshold,
            backgndThresRnd2=backgndThresRnd2,
            cellProportionAboveBackgnd=cellProportionAboveBackgnd,
            cellProportionBelowBackgndMode=cellProportionBelowBackgndMode,
            cellProportionForModeDetection=cellProportionForModeDetection,
            hashThreshold=hashThreshold,
            bgndQuantileThreshold=bgndQuantileThreshold,
            bgndQuantileThresRnd2=bgndQuantileThresRnd2,
            outputFormat=outputFormat)

        snacsObj=imputeMissingMutations(snacsObj,verbose=verbose)

        snacsObj=clusterCellsWithSNPdata(snacsObj,outputFormat="none")
            
        snacsObj=makeSnacsCall(snacsObj,
            minSampleSize=minSampleSize,
            cellProportionAboveBackgnd=cellProportionAboveBackgnd,
            minClustSize=minClustSize,
            clustComparePValue=clustComparePValue,
            maxClustSampleSize=maxClustSampleSize,
            clustCompareMethod=clustCompareMethod,
            maxClustSizeRnd2=maxClustSizeRnd2,
            dataTypeRnd2=dataTypeRnd2,
            cbsAlpha=cbsAlpha)
            
        for (snacsCallRnd in c("snacsRnd1","snacsRnd2")) {
            snacsCallThis=snacsObj$annCell[,snacsCallRnd]
            k=which(!snacsObj$annHash$hashNames%in%snacsCallThis)
            if (length(k)!=0) {
                if (bgndThresDetMethod=="two modes") {
                    warning(paste0('Cells cannot be assigned hash(es) ',paste0(snacsObj$annHash$hashNames[k],collapse=', '),'in ',snacsCallRnd,'. Try running function getBestSNPs with bgndThresDetMethod as "manual"'))
                } else {
                    done=F
                }
            } else {
                x=table(snacsCallThis[which(snacsCallThis%in%snacsObj$annHash$hashNames)])
                k=which(x<minSampleSize)
                if (length(k)!=0) {
                    if (bgndThresDetMethod=="two modes") {
                        warning(paste0('Not enough cells can be assigned hash(es) ',paste0(names(x)[k],collapse=', '),'in ',snacsCallRnd,'. Try running function getBestSNPs with bgndThresDetMethod as "manual" or by lowering parameter minSampleSize'))
                    } else {
                        done=F
                    }
                }
            }
        }
        if (done) break
    }

    clustObj=SNACS::createHeatmap(snacsObj,cell_anno_var=c("snacsRnd2","clustSplitRnd2","snacsRnd1","clustSplitRnd1","clustBestSNPs_hclust",snacsObj$annHash$hashNames),cell_anno_name=c("snacsRnd2","splitRnd2","snacsRnd1","splitRnd1","bestSNPsCluster",snacsObj$annHash$hashNames),col_dend=TRUE,row_dend=FALSE,outputFileName=paste0("heatmap_pvBestWithAnn_",snacsObj$exptName),outputFormat=outputFormat)

    invisible(snacsObj)
}

###########################################################
###########################################################
#' Select SNPs that best separate the hashes
#'
#' Group cells into constituent samples based on hash antibody data. Select the SNPs that best separate the groups.
#'
#' @param snacsObj SNACSList object
#' @param numSNP Integer. Number of pairs of SNPs (high in one sample / low in second sample and vice versa) that best separate sample-pairs. There will be s times n(n-1) number of SNPs where s is the number of SNPs and n is the number of samples. Default is 3
#' @param minSampleSize Integer. Minimum number of cells that is estimated by SNACS to be in a constituent sample. Gives a warning if below this number. Default is 100
#' @param bgndThresDetMethod Character. Method for detecting threshold of background antibody distribution. Default is "automatic"
#' @param backgndThreshold Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell for making first round of SNACS calls. Default is 0.95. Range is 0-1
#' @param backgndThresRnd2 Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell for making second round of SNACS calls. Default is 0.75. Range is 0-1
#' @param cellProportionAboveBackgnd Numeric. Proportion of cells in a cluster that have to have signal above the background of a hash, determined by "backgndThreshold" parameter, for the cluster to be assigned to the hash. Default is 0.5. Range is 0-1
#' @param cellProportionBelowBackgndMode Numeric. Maximum proportion of cells which can be below the mode of the estimated hash background distribution. Default is 0.6. Range is 0-1
#' @param cellProportionForModeDetection Numeric. Proportion of cells used to estimate mode of the background distribution. Used only if "cellProportionBelowBackgndMode" threshold is not met; otherwise, all cells are used. Default is 0.75. Range is 0-1
#' @param hashThreshold Numeric. Threshold of the hash value above which the antibody will be considered to be expressed in a cell. If NA then backgndThreshold is used. Default is 0.5
#' @param bgndQuantileThreshold Numeric. Threshold of the hash value above which the antibody will be considered to be expressed in a cell in first round of SNACS call. Default is NA. Used if bgndThresDetMethod = "manual"
#' @param bgndQuantileThresRnd2 Numeric. Threshold of the hash value above which the antibody will be considered to be expressed in a cell in second round of SNACS call. Default is NA. Used if bgndThresDetMethod = "manual"
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @param verbose Logical. Prints information when running the method. Default is FALSE
#' @return A SNACSList object
#' @export
getBestSNPs=function(snacsObj,numSNP=3,minSampleSize=100,bgndThresDetMethod=c("automatic","manual","default modes","two modes"),backgndThreshold=0.95,backgndThresRnd2=0.75,cellProportionAboveBackgnd=0.5,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,hashThreshold=0.5,bgndQuantileThreshold=NA,bgndQuantileThresRnd2=NA,outputFormat=c("","pdf","png"),verbose=FALSE) {
    
    bgndThresDetMethod=bgndThresDetMethod[1]
    
    ## -----------------------------------
    if (is.na(match("filtered",snacsObj[["processLevels"]]))) stop("Run filterData() before selecting best SNPs\n")

    ## -----------------------------------
    snacsObj=clusterSampleWithAntibodyData(snacsObj,minSampleSize=minSampleSize,bgndThresDetMethod=bgndThresDetMethod,backgndThreshold=backgndThreshold,cellProportionBelowBackgndMode=cellProportionBelowBackgndMode,cellProportionForModeDetection=cellProportionForModeDetection,hashThreshold=hashThreshold,bgndQuantileThreshold=bgndQuantileThreshold,bgndQuantileThresRnd2=bgndQuantileThresRnd2)
    generateAntibodyDensityPlot(snacsObj,backgndThreshold=backgndThreshold,cellProportionBelowBackgndMode=cellProportionBelowBackgndMode,cellProportionForModeDetection=cellProportionForModeDetection,outputFormat=outputFormat)

    ## Association of mutation status of each SNP with each sample cluster pair
    samUniq=snacsObj$annHash$hashNames
    n=nrow(snacsObj$annSNP)*length(samUniq)*(length(samUniq)-1)/2
    tmp=rep(NA,n); tmpC=rep("",n)
    stat=data.frame(id=tmpC,sam1=tmpC,sam2=tmpC,stat=tmp,pv=tmp,propMutIn1stSample=tmp)
    k=1
    optionThis=getOption("warn")
    options(warn=-1)
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
    options(warn=optionThis)
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
#' Cluster cells based on SNP data
#'
#' Cluster cells, where the number of clusters is the number of constituent samples, using the mutation data of select SNPs. Generate heatmap of the mutation data.
#'
#' @param snacsObj SNACSList object
#' @param outputFormat Character. Output file type. "" outputs to the standard output. Default is "none" which does not plot anything
#' @param verbose Logical. Prints information when running the method. Default is FALSE
#' @return A SNACSList object
#' @export
clusterCellsWithSNPdata=function(snacsObj,outputFormat=c("none","","pdf","png"),verbose=FALSE) {
    outputFormat=outputFormat[1]

    ## -----------------------------------
    if (is.na(match("bestSNPs",snacsObj[["processLevels"]]))) stop("Run getBestSNPs() before clustering SNPs\n")
    if (is.na(match("imputed",snacsObj[["processLevels"]]))) stop("Run imputeMissingMutations() before clustering SNPs\n")

    clustMethod=c("hclust","skmean")
    clustMethod=clustMethod[1]

    pvSnpThres=NA
    cellClusterFileName=NA
    subsetSnpFlag=""
    subsetCellFlag=""
    
    outputFileName=paste0("heatmap",ifelse(clustMethod=="skmean","_skmean",""),"_bestSNPs_",snacsObj$exptName)
    if (outputFormat=="") h_title="Heatmap of best SNPs" else h_title=NULL
    
    cell_anno_var_bestSNPs=snacsObj$annHash$hashNames
    
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
    names(annCellThis)[match(paste0("cluster_",clustMethod),names(annCellThis))]=paste0("clustBestSNPs_",clustMethod)

    snacsObj[["mut"]]=datThis
    snacsObj[["hashes"]]=hashesThis
    snacsObj[["annCell"]]=annCellThis
    snacsObj[["annSNP"]]=annSNPthis

    ################################################
    clustObj=createHeatmap(snacsObj,cell_anno_var=cell_anno_var_bestSNPs,col_dend=T,row_dend=F,h_title=h_title,outputFormat=outputFormat,outputFileName=outputFileName)
    snacsObj[["hclustObj_bestSNPs"]]=clustObj$colClust

    ################################################
    invisible(snacsObj)
}

####################################################################
####################################################################
#' Make SNACS calls
#'
#' The SNACS calls are made in two rounds. In round 1, the cells are clustered based on select SNPs. All cells in a cluster are assigned to a sample if > 50% of cells from that cluster have a hash expression that exceeds the 95th percentile of the background distribution for that sample. Clusters assigned to multiple hashes are designated as multiplets. In round 2, the multiplet detection is refined by segmenting the antibody data and assigning a segment to a sample when it has > 50% cells that exceeds the 75th percentile of the hash background. Only narrow segments (<=100 cells) are considered in this round. The results are stored in the columns "snacsRnd1" and "snacsRnd2" in the "annCell" table of the SNACS object
#'
#' @param snacsObj SNACSList object
#' @param minSampleSize Integer. Minimum number of cells that is estimated by SNACS to be in a constituent sample. Gives a warning if below this number. Default is 100
#' @param cellProportionAboveBackgnd Numeric. Proportion of cells in a cluster that have to have signal above the background of a hash, determined by "backgndThreshold" parameter, for the cluster to be assigned to the hash. Default is 0.5. Range is 0-1
#' @param minClustSize Numeric. Minimum number of cells required to be in a cluster. Default is 2
#' @param clustComparePValue Numeric. P-value threshold to compare cluster pairs. Default is 10^-5
#' @param maxClustSampleSize Numeric. Maximum number of cells in a cluster that will be used to compare. If more, then this number of cells will be sampled. Default is Inf
#' @param clustCompareMethod Character. Test used to compare clusters. Default is t-test
#' @param maxClustSizeRnd2 Integer. Maximum number of cells required to be in a cluster for making second round of SNACS calls. Default is 100
#' @param dataTypeRnd2 Character. Type of data to be used for splitting the cells when making second round of SNACS calls. Default is euclidean distance of the cells to the cluster means from first round of SNACS calls
#' @param cbsAlpha Numeric. Significance level passed to DNAcopy::segment to for splitting the cells when making second round of SNACS calls. Default is 0.1
#' @param verbose Logical. Prints information when running the method. Default is FALSE
#' @return A SNACSList object
#' @export
makeSnacsCall=function(snacsObj,minSampleSize=100,cellProportionAboveBackgnd=0.5,minClustSize=2,clustComparePValue=10^-5,maxClustSampleSize=Inf,clustCompareMethod=c("t","hotelling"),maxClustSizeRnd2=100,dataTypeRnd2=c("euclidean","sum of squares","log2 euclidean","log2 sum of squares"),cbsAlpha=0.1,verbose=FALSE) {
    clustCompareMethod=clustCompareMethod[1]
    dataTypeRnd2=dataTypeRnd2[1]

    #Make second round of SNACS calls to detect narrow multiplet regions. Always TRUE
    makeSnacsCallRnd2=T

    ## -----------------------------------
    if (is.na(match("hclustObj_bestSNPs",names(snacsObj)))) stop("Run clusterCellsWithSNPdata() before making SNACS calls\n")

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
    
    ## --------------------------------------
    ## --------------------------------------
    
    annHashBackgnd=snacsObj$annHashBackgnd

    ###########################################################
    ## Make SNACS Round 1 calls
    
    clustInfo$snacsRnd1="Multiplet"
    grp=clustInfo$clustId
    grpUniq=sort(unique(grp))
    for (gId in 1:length(grpUniq)) {
        j=which(grp==grpUniq[gId])
        x=apply(hashMat[j,],2,mean,na.rm=T)
        clustInfo$snacsRnd1[j]=names(x)[which.max(x)]
    }

    snacsRnd1=clustInfo$snacsRnd1
    
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
            if (round(mean(hashMatThis[j,sampleId]>annHashBackgnd$thres[which(annHashBackgnd$hashNames==sampleId)],na.rm=T),2)>=cellProportionAboveBackgnd) {inSamples=c(inSamples,sampleId)}
        }
        snacsRnd1[j]=paste(inSamples,collapse="_")
    }
    
    ## Mark cluster splits
    j=clustObjThis$order
    grp=as.integer(as.factor(x[j]))
    j=which(diff(grp)!=0)
    jj=j
    j1=c(1,j+1); j2=c(j,length(grp))
    if (F) {
        for (j in 1:length(j1)) {
            if (any(grp[j1[j]:j2[j]]!=grp[j1[j]])) print(j)
        }
    }
    clustSplitRnd1=rep(1,length(grp))
    for (j in 1:length(j1)) {
        if (j%%2==0) clustSplitRnd1[j1[j]:j2[j]]=2
    }
    clustSplitRnd1=clustSplitRnd1[match(clustInfo$id,clustObjThis$label[clustObjThis$order])]
    subClustRnd1=x

    ## Distance to centroid
    grp=snacsRnd1
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
    ## Make SNACS Round 2 calls
    
    if (makeSnacsCallRnd2) {
        ## Detect and assign narrow regions as multiplets
        ## Make splits based on cbs on distance to centroids
        
        set.seed(12345)
        cellId=which(snacsRnd1%in%snacsObj$annHash$hashNames)
        #x1=hashMatThis[clustObjThis$order,]
        x1=t(dist2centroidMat[,clustObjThis$order])
        x2=DNAcopy::CNA(genomdat=x1,
        chrom=rep(1,nrow(x1)),maploc=1:nrow(x1),
                          data.type="logratio",sampleid=colnames(x1))
        x3=DNAcopy::segment(x2,alpha=cbsAlpha,verbose=0)
        x3=x3$output
        x=rep(NA,nrow(x1))
        names(x)=clustObjThis$label[clustObjThis$order]
        grp=snacsRnd1[clustObjThis$order]
        grpUniq=snacsObj$annHash$hashNames
        for (gId in 1:length(grpUniq)) {
            j=which(grp==grpUniq[gId])
            if (length(j)!=0) {
                j1=range(j,na.rm=T)
                x2=rep(NA,nrow(x1))
                p1=1
                for (k in which(x3$ID==grpUniq[gId])) {
                    x2[x3$loc.start[k]:x3$loc.end[k]]=p1
                    p1=p1+1
                }
                x[j1[1]:j1[2]]=paste0(grpUniq[gId],"_",x2[j1[1]:j1[2]])
            }
        }
        x=x[match(clustObjThis$label,names(x))]
        x=x[cellId]

        grp=snacsRnd1
        for (k in unique(x[!is.na(x)])) {
            j=cellId[which(x==k)]
            if (length(j)>maxClustSizeRnd2) next
            inSamples=c()
            for (sampleId in snacsObj$annHash$hashNames) {
                if (round(mean(hashMatThis[j,sampleId]>annHashBackgnd$thresRnd2[which(annHashBackgnd$hashNames==sampleId)]),2)>=cellProportionAboveBackgnd) {inSamples=c(inSamples,sampleId)}
            }
            grp[j]=paste(inSamples,collapse="_")
        }
        snacsRnd2=grp

        ## Mark cluster splits
        cellId=1:nrow(clustInfoThis)
        cellId=which(snacsRnd1%in%snacsObj$annHash$hashNames)
        z1=rep("0",length(snacsRnd1)); z1[cellId]=x
        j=clustObjThis$order
        grp=as.integer(as.factor(z1[j]))
        j=which(diff(grp)!=0)
        jj=j
        j1=c(1,j+1); j2=c(j,length(grp))
        if (F) {
            for (j in 1:length(j1)) {
                if (any(grp[j1[j]:j2[j]]!=grp[j1[j]])) print(j)
            }
        }
        clustSplitRnd2=rep(1,length(grp))
        for (j in 1:length(j1)) {
            if (j%%2==0) clustSplitRnd2[j1[j]:j2[j]]=2
        }
        clustSplitRnd2=clustSplitRnd2[match(clustInfo$id,clustObjThis$label[clustObjThis$order])]
        subClustRnd2=z1
    } ## end of makeSnacsCallRnd2

    ###########################################################
    ###########################################################
    
    snacsObj$annCell$snacsRnd1=""
    j=match(snacsObj$hclustObj_bestSNPs$label,snacsObj$annCell$id); j1=which(!is.na(j)); j2=j[j1]
    snacsObj$annCell$snacsRnd1[j2]=snacsRnd1[j1]
    snacsObj$annCell$snacsRnd1[snacsObj$snacsRnd1==""]=NA
    snacsObj$annCell$clustSplitRnd1=clustSplitRnd1
    snacsObj$annCell$subClustRnd1=subClustRnd1
    snacsObj$centroid=centroid
    snacsObj$dist2centroidMat=dist2centroidMat
    if (makeSnacsCallRnd2) {
        snacsObj$annCell$snacsRnd2=""
        snacsObj$annCell$snacsRnd2[j2]=snacsRnd2[j1]
        snacsObj$annCell$snacsRnd2[snacsObj$snacsRnd2==""]=NA
        snacsObj$annCell$clustSplitRnd2=clustSplitRnd2
        snacsObj$annCell$subClustRnd2=subClustRnd2
    }
    
    if (verbose) cat('"snacs" column(s) added to "annCell" table in SNACS object\n\n',sep="")

    if (makeSnacsCallRnd2) {snacsCall=snacsRnd2
    } else {snacsCall=snacsRnd1}
    if (verbose) print(table(snacsCall,cluster=clustInfo[match(snacsObj$hclustObj_bestSNPs$label,clustInfo$id),"clustId"],exclude=NULL,dnn=list("SNACS calls",paste0("Cell clusters with best SNPs"))))
    rm(snacsCall)
    
    if (verbose) {
        for (snacsCallRnd in c("snacsRnd1","snacsRnd2")) {
            snacsCallThis=snacsObj$annCell[,snacsCallRnd]
            k=which(!snacsObj$annHash$hashNames%in%snacsCallThis)
            if (length(k)!=0) {warning(paste0('Cells cannot be assigned hash(es) ',paste0(snacsObj$annHash$hashNames[k],collapse=', '),'in ',snacsCallRnd,'. Try running function getBestSNPs with bgndThresDetMethod as "two modes" or "manual"'))
            } else {
                x=table(snacsCallThis[which(snacsCallThis%in%snacsObj$annHash$hashNames)])
                k=which(x<minSampleSize)
                if (length(k)!=0) {warning(paste0('Not enough cells can be assigned hash(es) ',paste0(names(x)[k],collapse=', '),'in ',snacsCallRnd,'. Try running function getBestSNPs with bgndThresDetMethod as "two modes" or "manual" or by lowering parameter minSampleSize'))}
            }
        }
    }

    invisible(snacsObj)
}

####################################################################
####################################################################
#' Group cells based on hash antibody data
#'
#' A cell is assigned to a specific sample if the antibody expression of that sample is expressed highly for that cell. These initial classifications are used to select SNPs for making SNACS calls.
#'
#' @param snacsObj SNACSList object
#' @param minSampleSize Integer. Minimum number of cells that is estimated by SNACS to be in a constituent sample. Gives a warning if below this number. Default is 100
#' @param bgndThresDetMethod Character. Method for detecting threshold of background antibody distribution. Default is "automatic"
#' @param backgndThreshold Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell for making first round of SNACS calls. Default is 0.95. Range is 0-1
#' @param backgndThresRnd2 Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell for making second round of SNACS calls. Default is 0.75. Range is 0-1
#' @param cellProportionBelowBackgndMode Numeric. Maximum proportion of cells which can be below the mode of the estimated hash background distribution. Default is 0.6. Range is 0-1
#' @param cellProportionForModeDetection Numeric. Proportion of cells used to estimate mode of the background distribution. Used only if "cellProportionBelowBackgndMode" threshold is not met; otherwise, all cells are used. Default is 0.75. Range is 0-1
#' @param hashThreshold Numeric. Threshold of the hash value above which the antibody will be considered to be expressed in a cell. If NA then backgndThreshold is used. Default is 0.5
#' @param bgndQuantileThreshold Numeric. Threshold of the hash value above which the antibody will be considered to be expressed in a cell in first round of SNACS call. Default is NA. Used if bgndThresDetMethod = "manual"
#' @param bgndQuantileThresRnd2 Numeric. Threshold of the hash value above which the antibody will be considered to be expressed in a cell in second round of SNACS call. Default is NA. Used if bgndThresDetMethod = "manual"
#' @param verbose Logical. Prints information when running the method. Default is FALSE
#' @export
clusterSampleWithAntibodyData=function(snacsObj,minSampleSize=100,bgndThresDetMethod=c("automatic","manual","default modes","two modes"),backgndThreshold=0.95,backgndThresRnd2=0.75,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,hashThreshold=0.5,bgndQuantileThreshold=NA,bgndQuantileThresRnd2=NA,verbose=FALSE) {
    bgndThresDetMethod=bgndThresDetMethod[1]
    
    
    hashMat=t(snacsObj$hashes)

    methodVec=bgndThresDetMethod; if (bgndThresDetMethod=="automatic") {methodVec=c("default modes","two modes")}
    tmp=rep(NA,length=nrow(snacsObj$annHash))
    tmpC=vector(mode="character",length=nrow(snacsObj$annHash))
    annHashBackgnd=data.frame(hashNames=snacsObj$annHash$hashNames,thres=tmp,thresRnd2=tmp,mode=tmp,adjust=tmp,bw=tmpC,method=tmpC)
    for (k in 1:ncol(hashMat)) {
        if (bgndThresDetMethod=="manual") {
            annHashBackgnd[k,c("thres","thresRnd2")]=c(bgndQuantileThreshold,bgndQuantileThresRnd2)
            annHashBackgnd$bw[k]=NA
            annHashBackgnd$method[k]=bgndThresDetMethod
        } else {
            for (methodThis in methodVec) {
                annHashBackgnd$method[k]=methodThis
                if (methodThis=="default modes") {
                    annHashBackgnd$bw[k]="SJ"
                    x=stats::density(hashMat[,k],bw=annHashBackgnd$bw[k],adjust=1,na.rm=T)
                    k1=which.max(x$y)
                    if (stats::quantile(hashMat[,k],probs=cellProportionBelowBackgndMode)<x$x[k1]) k1=which.max(x$y[1:stats::quantile(1:(k1-1),probs=cellProportionForModeDetection)])
                    annHashBackgnd[k,c("mode","adjust")]=c(x$x[k1],1)
                } else if (methodThis=="two modes") {
                    annHashBackgnd$bw[k]="nrd0"
                    annHashBackgnd[k,c("mode","adjust")]=getHashModeInfo(mm=hashMat[,k],bw=annHashBackgnd$bw[k])
                    x=stats::density(hashMat[,k],bw=annHashBackgnd$bw[k],adjust=annHashBackgnd$adjust[k],na.rm=T)
                }

                xx=annHashBackgnd$mode[k]
                x2=x$x[which(x$x<xx)]-xx; x2=c(x2,-x2); x2=x2+xx
                hashBackgndData=x2
                hashBackgndECDF=stats::ecdf(hashBackgndData)
                x2=hashMat[which(hashMat[,k]<xx),k]-xx; x2=c(x2,-x2); x2=x2+xx
                ecdfbg=1-hashBackgndECDF(x2)
                annHashBackgnd[k,c("thres","thresRnd2")]=c(stats::quantile(x2,probs=c(backgndThreshold,backgndThresRnd2)))
                
                if (sum(hashMat[,k]>annHashBackgnd$thres[k],na.rm=T)>=minSampleSize) break
            }
        }
        
    }
    
    snacsCallThis=rep("",nrow(snacsObj$annCell))
    j=1:length(snacsCallThis)
    for (k in 1:nrow(snacsObj$annHash)) {
        jj=which(hashMat[,k]>annHashBackgnd$thres[k])
        if (length(jj)!=0) snacsCallThis[jj]=paste0(snacsCallThis[jj],"_",annHashBackgnd$hashNames[k])
    }
    snacsCallThis=sub("^_", "",snacsCallThis)
    k=which(!snacsObj$annHash$hashNames%in%snacsCallThis)
    if (length(k)!=0) {warning(paste0('Cells cannot be assigned hash(es) ',paste0(snacsObj$annHash$hashNames[k],collapse=', '),'. Try running function getBestSNPs with bgndThresDetMethod as "two modes" or "manual"'))
    } else {
        x=table(snacsCallThis[which(snacsCallThis%in%snacsObj$annHash$hashNames)])
        k=which(x<minSampleSize)
        if (length(k)!=0) {warning(paste0('Not enough cells can be assigned hash(es) ',paste0(names(x)[k],collapse=', '),'. Try running function getBestSNPs with bgndThresDetMethod as "two modes" or "manual" or by lowering parameter minSampleSize'))}
    }

    snacsCallThis=rep("",nrow(snacsObj$annCell))
    j=1:length(snacsCallThis)
    for (k in 1:nrow(snacsObj$annHash)) {
        if (!is.na(hashThreshold)) {
            jj=which(hashMat[,k]>=hashThreshold)
        } else {
            jj=which(hashMat[,k]>annHashBackgnd$thres[k])
        }
        if (length(jj)!=0) snacsCallThis[jj]=paste0(snacsCallThis[jj],"_",annHashBackgnd$hashNames[k])
    }
    snacsCallThis=sub("^_", "",snacsCallThis)
    
    snacsObj$annHashBackgnd=annHashBackgnd
    snacsObj$annCell$hashClust=snacsCallThis

    invisible(snacsObj)
}

####################################################################
####################################################################
#' Generate hash antibody background density plots
#'
#' The density plots are helpful is understanding the background distribution of a hash antibody. The plots are saved in "../output" folder of a "pdf" or "png" format is specified. There are 3 figures generated for each sample. The top figure shows the distribution of a hash antibody. The background distribution, shown in red, is estimated by generating a symmetric distribution by reflecting the data to the left of the left mode about that mode. The green lines mark the median and 95th percentile of the background distribution. The middle figure shows just the estimated background distribution. The bottom figure is the histogram of the hash antibody data.
#'
#' @param snacsObj SNACSList object
#' @param backgndThreshold Numeric. Threshold of the background antibody distribution of a hash above which the antibody will be considered to be expressed in a cell for making first round of SNACS calls. Default is 0.95. Range is 0-1
#' @param cellProportionBelowBackgndMode Numeric. Maximum proportion of cells which can be below the mode of the estimated hash background distribution. Default is 0.6. Range is 0-1
#' @param cellProportionForModeDetection Numeric. Proportion of cells used to estimate mode of the background distribution. Used only if "cellProportionBelowBackgndMode" threshold is not met; otherwise, all cells are used. Default is 0.75. Range is 0-1
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @export
generateAntibodyDensityPlot=function(snacsObj,backgndThreshold=0.95,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,outputFormat=c("","pdf","png")) {
    #snacsObj=snacsObj; backgndThreshold=0.95; cellProportionBelowBackgndMode=0.6; cellProportionForModeDetection=0.75; outputFormat=""
    
    outputFormat=outputFormat[1]

    ## -----------------------------------
    if (outputFormat=="") {
        plotInfo=list(cexMain=1,cexLab=1,cexAxis=1)
    } else {
        dirResult="../output/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
        dirResult="../output/hashDensityPlot/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
        #dirResult="../output/hashDensityPlot/"
        plotInfo=list(cexMain=2,cexLab=1.5,cexAxis=1.5)
    }
    
    hashMat=t(snacsObj$hashes)
    ## --------------------------------------

    annHashBackgnd=snacsObj$annHashBackgnd
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
        
        bwMethod=annHashBackgnd$bw[k]
        xx=annHashBackgnd$mode[k]

        x=stats::density(hashMat[,k],bw=bwMethod,adjust=annHashBackgnd$adjust[k],na.rm=T)
        xlim=range(x$x); ylim=range(x$y); ylim=NULL; ylim=c(0,1)
        xlim=max(abs(range(x$x))); xlim=c(-xlim,xlim)
        x2=x$x[which(x$x<xx)]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBGSymData=x2
        hashBackgndECDF=stats::ecdf(hashBGSymData)
        x2=hashMat[which(hashMat[,k]<xx),k]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBGData=x2
        ecdfbg=1-hashBackgndECDF(hashBGData)
        vertLine=c(stats::median(x2),stats::quantile(x2,probs=backgndThreshold)); vertLineLab=c("med",backgndThreshold)
        graphics::plot(x,xlim=xlim,ylim=ylim,main=paste0(snacsObj$exptName,": ",colnames(hashMat)[k]),xlab="Hash",cex.main=plotInfo$cexMain,cex.lab=plotInfo$cexLab,cex.axis=plotInfo$cexAxis)
        graphics::lines(stats::density(x2,bw=bwMethod,adjust=annHashBackgnd$adjust[k],na.rm=T),col="red"); graphics::abline(v=vertLine,col="green")
        graphics::plot(stats::density(x2,bw=bwMethod,adjust=annHashBackgnd$adjust[k],na.rm=T),xlim=xlim,ylim=ylim,main="",xlab="Hash background",col="red",cex.main=plotInfo$cexMain,cex.lab=plotInfo$cexLab,cex.axis=plotInfo$cexAxis); graphics::abline(v=vertLine,col="green")
        graphics::axis(side=3,at=vertLine,labels=vertLineLab,cex.axis=plotInfo$cexAxis,las=3,col="green")
        x=hashMat[,k]
        graphics::hist(x,freq=T,breaks=diff(range(x))/.1,xlim=xlim,ylim=ylimHist,main=paste0(snacsObj$exptName,": ",colnames(hashMat)[k]),xlab="Hash",ylab="Count",cex.main=plotInfo$cexMain,cex.lab=plotInfo$cexLab,cex.axis=plotInfo$cexAxis)

        if (outputFormat!="") grDevices::dev.off()
    }
}

####################################################################
####################################################################
#' Detect density function modes
#'
#' Find the modes of the density function of the hash antibody distribution for background detection.
#'
#' @param x Numeric. Hash antibody data
#' @param bw Character. Smoothing bandwidth parameter passed to density function
#' @param adjust Character. Bandwidth adjustment factor parameter passed to density function
#' @param signifi Numeric. Number of significant digits of the density modes, Used to order the modes
#' @param from Numeric. left-most point value passed to density function
#' @param to Numeric. right-most point value passed to density function
#' @return A vector of modes
#' @export
getModes=function(x,bw,adjust,signifi,from,to) {
    den=stats::density(x, kernel=c("gaussian"),bw=bw,adjust=adjust,from=from,to=to)
    den.s=stats::smooth.spline(den$x, den$y, all.knots=TRUE, spar=0.1)
    s.1=stats::predict(den.s, den.s$x, deriv=1)
    s.0=stats::predict(den.s, den.s$x, deriv=0)
    den.sign=sign(s.1$y)
    a=c(1,1+which(diff(den.sign)!=0))
    b=rle(den.sign)$values
    df=data.frame(a,b)
    df = df[which(df$b %in% -1),]
    modes=s.1$x[df$a]
    dens=s.0$y[df$a]
    df2=data.frame(modes=modes,density=dens)
    df2$sig=signif(df2$density,signifi)
    df2=df2[with(df2,order(-sig)),]
    
    df2[order(df2$modes),]
}

####################################################################
####################################################################
#' Get mode information of density function
#'
#' Get mode information of hash antibody distribution for background detection.
#'
#' @param mm Numeric. Hash antibody data
#' @param bw Character. Smoothing bandwidth parameter passed to density function
#' @param verbose Logical. Prints information when running the method. Default is FALSE
#' @return Mode and adjustment values
#' @export
getHashModeInfo=function(mm,bw="nrd0",verbose=FALSE) {
    adjust=adjustThis=adjustPrev=1; offset=0.25
    adjustMax=diff(range(mm))/2
    offsetMM=0.001*diff(range(mm))

    dirn=""
    while (T) {
        modeInfo=getModes(mm,bw=bw,adjust=adjust,2,min(mm)-offsetMM,max(mm)+offsetMM)

        if (nrow(modeInfo)==2) {
            break
        } else if (nrow(modeInfo)>2) {
            adjustThis=adjustThis+offset
            if (adjustThis==1 | adjustThis>adjustMax | dirn=="dn") {
                #if (dirn=="dn") cat("dirn is down\n")
                if (verbose) cat('Getting more than two modes for hash density function. Check hash density plots and SNACS heatmaps. Might try running function getBestSNPs with bgndThresDetMethod="manual"\n')
                adjustThis=adjust=adjustPrev
                modeInfo=getModes(mm,bw=bw,adjust=adjust,2,min(mm)-offsetMM,max(mm)+offsetMM)
                break
            }
            dirn="up"
        } else {
            adjustThis=adjustThis-offset
            if (adjustThis==1 | adjustThis<=0 | dirn=="up") {
                #if (dirn=="up") cat("dirn is up\n")
                if (verbose) cat('Getting single mode for hash density function. Check hash density plots and SNACS heatmaps. Might try running function getBestSNPs with bgndThresDetMethod="manual"\n')
                adjustThis=adjust=adjustPrev
                modeInfo=getModes(mm,bw=bw,adjust=adjust,2,min(mm)-offsetMM,max(mm)+offsetMM)
                break
            }
            dirn="dn"
        }
        adjust=adjustThis
    }
    
    modeInfo=modeInfo[order(modeInfo$density,decreasing=T),]
    modeInfo=modeInfo[1:min(2,nrow(modeInfo)),]
    modeInfo=modeInfo[order(modeInfo$mode),]
    
    out=c(modeInfo$modes[1],adjust)
    names(out)=c("mode","adjust")
    
    out
}

###########################################################
###########################################################
