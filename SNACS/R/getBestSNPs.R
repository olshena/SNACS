####################################################################
####################################################################
#' Select SNPs that best separate the hashes
#'
#' Cluster cells, where the number of clusters is the number of hashes, using the mutation data. Rank the SNPs in the order that best separates into those clusters. Then select the best SNPs that separate the cells into those clusters and re-cluster the cells using those SNPs. Generate heatmaps of the mutation data with the ranked and best SNPs.
#'
#' @param snacsObj SNACSList object
#' @param cell_anno_var Character vector. Cell variables for which color bars are to be added to the heatmap. Default is NULL
#' @param clustMethodForRankedSNPs Character. Clustering method for ranking SNPs. Options are "hclust" and "skmean". Default is "hclust"
#' @param clustMethodForBestSNPs Character. Clustering method for selecting best SNPs. Options are "hclust" and "skmean". Default is "hclust"
#' @param pvSnpThres Numeric. P-value cutoff to select best SNPs. Default is 10^-320
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @return A SNACSList object
#' @export
getBestSNPs=function(snacsObj,cell_anno_var=NULL,clustMethodForRankedSNPs=c("hclust","skmean"),clustMethodForBestSNPs=c("hclust","skmean"),pvSnpThres=10^-320,outputFormat=c("","pdf","png")) {
    if (is.na(match("imputed",snacsObj[["processLevels"]]))) stop("Run imputeMissingMutations() before selecting SNPs for making hash calls\n")
    if (is.null(cell_anno_var)) cell_anno_var_rankedSNPs=snacsObj$annHash$hashNames
    cell_anno_var_bestSNPs=c(cell_anno_var_rankedSNPs,paste0("clustRankedSNPs_",clustMethodForRankedSNPs[1]))
    snacsObj=getRankedSNPs.internal(snacsObj,cell_anno_var=cell_anno_var_rankedSNPs,clustMethod=clustMethodForRankedSNPs,outputFormat=outputFormat)
    snacsObj=getBestSNPs.internal(snacsObj,cell_anno_var=cell_anno_var_bestSNPs,clustMethod=clustMethodForBestSNPs,pvSnpThres=pvSnpThres,outputFormat=outputFormat)
    invisible(snacsObj)
}

###########################################################
###########################################################
#' Rank SNPs in the order that best separates the hashes
#'
#' Cluster cells, where the number of clusters is the number of hashes, using the mutation data. Then rank the SNPs in the order that best separates into those clusters. Generate heatmap of the mutation data with the ranked SNPs.
#'
#' @param snacsObj SNACSList object
#' @param cell_anno_var Character vector. Cell variables for which color bars are to be added to the heatmap. Default is NULL
#' @param clustMethod Character. Clustering method for ranking SNPs. Options are "hclust" and "skmean". Default is "hclust"
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @return A SNACSList object
#' @export
getRankedSNPs.internal=function(snacsObj,cell_anno_var,clustMethod=c("hclust","skmean"),outputFormat=c("","pdf","png")) {
    clustMethod=clustMethod[1]
    outputFormat=outputFormat[1]

    pvSnpThres=NA
    cellClusterFileName=NA
    subsetSnpFlag=""
    subsetCellFlag=""
    
    outputFileName=paste0("heatmap",ifelse(clustMethod=="skmean","_skmean",""),"_pvRanked_",snacsObj$exptName)
    if (outputFormat=="") h_title="Heatmap of ranked SNPs" else h_title=NULL
    
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
    
    j=apply(datThis,2,function(x) {mean(x==1)}); j=which(j>0 & j<1) ## Exclude cells with no or all mutations
    datThis=datThis[,j]
    hashesThis=hashesThis[,j]
    annCellThis=annCellThis[j,]
    
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
    }

    ################################################
    if (subsetSnpFlag!="") {
        if (subsetSnpFlag=="_knownHashSnpOrder") {
            tbl=utils::read.table(paste0("../heatmap/_knownHash/clustInfoSnp_knownHash_",snacsObj$exptName,".txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            i=match(tbl$id,annSNPthis$id)
            i=rev(i)
            datThis=datThis[i,]
            annSNPthis=annSNPthis[i,]
        }
    }
    
    ################################################
    j=apply(datThis,2,function(x) {mean(x==1)}); j=which(j>0 & j<1) ## Exclude cells with no or all mutations
    datThis=datThis[,j]
    hashesThis=hashesThis[,j]
    annCellThis=annCellThis[j,]

    ################################################
    names(annCellThis)[match(paste0("cluster_",clustMethod),names(annCellThis))]=paste0("clustRankedSNPs_",clustMethod)

    snacsObj[["mut"]]=datThis
    snacsObj[["hashes"]]=hashesThis
    snacsObj[["annCell"]]=annCellThis
    snacsObj[["annSNP"]]=annSNPthis

    ################################################
    createHeatmap(snacsObj,cell_anno_var,col_dend=T,row_dend=F,h_title=h_title,outputFormat=outputFormat,outputFileName=outputFileName)
    
    ################################################
    invisible(snacsObj)
}

###########################################################
###########################################################
#' Select SNPs that best separates the hashes
#'
#' Select the best SNPs that separate the cells into clusters, where the number of clusters is the number of hashes, using the mutation data. Then re-cluster the cells using the best SNPs. Generate heatmap of the mutation data with the best SNPs.
#'
#' @param snacsObj SNACSList object
#' @param cell_anno_var Character vector. Cell variables for which color bars are to be added to the heatmap. Default is NULL
#' @param clustMethod Character. Clustering method for selecting best SNPs. Options are "hclust" and "skmean". Default is "hclust"
#' @param pvSnpThres Numeric. P-value cutoff to select best SNPs. Default is 10^-320
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @return A SNACSList object
#' @export
getBestSNPs.internal=function(snacsObj,cell_anno_var,clustMethod=c("hclust","skmean"),pvSnpThres=10^-320,outputFormat=c("","pdf","png")) {
    clustMethod=clustMethod[1]
    cellClusterFileName=NA
    subsetSnpFlag=""
    subsetCellFlag=""
    outputFormat=outputFormat[1]
    
    outputFileName=paste0("heatmap",ifelse(clustMethod=="skmean","_skmean",""),"_pvBest_",snacsObj$exptName)
    if (outputFormat=="") h_title="Heatmap of best SNPs" else h_title=NULL

    ################################################
    ## Compute cosine cross-distances between the rows of matrices
    getSKmeansDist=function(x) {stats::as.dist(skmeans::skmeans_xdist(x))}
    distfun=getSKmeansDist
    linkMethod="ward.D2"

    ################################################
    datThis=snacsObj$mut
    hashesThis=snacsObj$hashes
    annCellThis=snacsObj$annCell
    annSNPthis=snacsObj$annSNP
    
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
        annSNPthis$minPValue_2Vs1Cluster=apply(pvMat,1,min,na.rm=T)
        annSNPthis$maxLOR_1VsOtherClusters=apply(lorMat,1,max,na.rm=T)
        i=rev(order(annSNPthis$minPValue_1VsOtherClusters))
    } else {
        if (F) {
            if (outputFormat!="") {
                dirResult="../output/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
                dirResult="../output/cellClusterPairPvalueHistogram/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
                grDevices::png(paste0(dirResult,"histogram_cellClusterPair_withRankedSNPs_",snacsObj$exptName,".png"),width=480,height=480)
            }
            x=annSNPthis$minPValue_1VsOtherClusters
            graphics::hist(x,freq=T,breaks=length(x)/100,main=paste0(snacsObj$exptName),xlab=paste0("1 cell cluster vs. other ",nrow(snacsObj$annHash)-1,": Minimum p-value"),ylab="Count",cex.main=2,cex.lab=1.5,cex.axis=1.5)
            if (outputFormat!="") grDevices::dev.off()
        }

        i=which(annSNPthis$minPValue_1VsOtherClusters<=pvSnpThres)
        if (length(i)==0) {
            cat("\n\nNo. of SNPs: ",nrow(annSNPthis),"\n",sep="")
            cat("SNP p-values:\n")
            print(signif(summary(annSNPthis$minPValue_1VsOtherClusters),2))
            stop(paste0("No SNPs with p-value <= ",pvSnpThres))
        }
    }
    datThis=datThis[i,]
    annSNPthis=annSNPthis[i,]
    
    j=apply(datThis,2,function(x) {mean(x==1)}); j=which(j>0 & j<1) ## Exclude cells with no or all mutations
    datThis=datThis[,j]
    hashesThis=hashesThis[,j]
    annCellThis=annCellThis[j,]
    
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
    }

    ################################################
    if (subsetSnpFlag!="") {
        if (subsetSnpFlag=="_knownHashSnpOrder") {
            tbl=utils::read.table(paste0("../heatmap/_knownHash/clustInfoSnp_knownHash_",snacsObj$exptName,".txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            i=match(tbl$id,annSNPthis$id)
            i=rev(i)
            datThis=datThis[i,]
            annSNPthis=annSNPthis[i,]
        }
    }

    ################################################
    i=match(annSNPthis$id,snacsObj$annSNP$id)
    snacsObj$mut=snacsObj$mut[i,]
    snacsObj$annSNP=snacsObj$annSNP[i,]
    
    ################################################
    j=apply(snacsObj$mut,2,function(x) {mean(x==1)}); j=which(j>0 & j<1) ## Exclude cells with no or all mutations
    snacsObj$mut=snacsObj$mut[,j]
    snacsObj$hashes=snacsObj$hashes[,j]
    snacsObj$annCell=snacsObj$annCell[j,]

    ################################################
    clustObj=createHeatmap(snacsObj,cell_anno_var,col_dend=T,row_dend=F,h_title=h_title,outputFormat=outputFormat,outputFileName=outputFileName)
    snacsObj[["hclustObj_bestSNPs"]]=clustObj$colClust
    
    ################################################
    invisible(snacsObj)
}

###########################################################
###########################################################
