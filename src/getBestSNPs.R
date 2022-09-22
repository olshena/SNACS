getBestSNPs=function(snacsObj,cell_anno_var=NULL,clustMethodForRankedSNPs=c("hclust","skmean"),clustMethodForBestSNPs=c("hclust","skmean"),pvSnpThres=0) {
    if (is.null(cell_anno_var)) cell_anno_var_rankedSNPs=snacsObj$annHash$hashNames
    cell_anno_var_bestSNPs=c(cell_anno_var_rankedSNPs,paste0("clustRankedSNPs_",clustMethodForRankedSNPs[1]))
    snacsObj=getRankedSNPs.internal(snacsObj,cell_anno_var=cell_anno_var_rankedSNPs,clustMethod=clustMethodForRankedSNPs)
    snacsObj=getBestSNPs.internal(snacsObj,cell_anno_var=cell_anno_var_bestSNPs,clustMethod=clustMethodForBestSNPs,pvSnpThres=pvSnpThres)
    invisible(snacsObj)
}

###########################################################
###########################################################
getRankedSNPs.internal=function(snacsObj,cell_anno_var,clustMethod=c("hclust","skmean")) {
    
    clustMethod=clustMethod[1]

    pvSnpThres=NA
    cellClusterFileName=NA
    subsetSnpFlag=""
    #subsetCellFlag="_allCells"
    subsetCellFlag=""
    
    outputfileName=paste0("heatmap",ifelse(clustMethod=="skmean","_skmean",""),"_pvRanked_",snacsObj$exptName)
    
    ################################################
    library(skmeans)
    library(RColorBrewer)
    library(marray)
    dirSrc="../src/"
    source(paste0(dirSrc,"heatmap4.R"))
    source(paste0(dirSrc,"heatmapRelated.R"))
    source(paste0(dirSrc,"generate_heatmap.R"))

    library(skmeans)
    getSKmeansDist=function(x) {as.dist(skmeans_xdist(x))}
    distfun=getSKmeansDist
    linkMethod="ward.D2"

    ## Distance functions available in "heatmapRelated.R"
    #getDist=getCosineDist
    #getDist=dist

    ################################################
    datThis=snacsObj$mut
    annSNPthis=snacsObj$annSNP
    annCellThis=snacsObj$annCell
    hashesThis=snacsObj$hashes
    if (F) {
        annCellThis$hash[annCellThis$hash==""]=NA
        annCellThis$hashKnown=annCellThis$hash
        annCellThis$hashKnown[!annCellThis$hash%in%c("CD45-26","CD45-27","CD45-28")]=NA
        annCellThis$hashUnknown=annCellThis$hash
        annCellThis$hashUnknown[annCellThis$hash%in%c("CD45-26","CD45-27","CD45-28")]=NA
    }
    
    ################################################
    if (T) {
        if (is.na(pvSnpThres)) {
            if (clustMethod=="skmean") {
                library(skmeans)
                library(cluster)
                clustCell=skmeans(t(datThis),k=nrow(snacsObj$annHash))
                silCell=silhouette(clustCell)
                clustCell=clustCell$cluster[match(annCellThis$id,names(clustCell$cluster))]
                silCell=silCell[match(annCellThis$id,rownames(silCell)),]
                annCellThis$cluster_skmean=clustCell
                annCellThis$cluster_silhouette=silCell[,"sil_width"]
                grp=annCellThis$cluster_skmean
            } else if (!is.na(cellClusterFileName)) {
                clustCell=read.table(cellClusterFileName, sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=-1)
                j=match(annCellThis$id,clustCell$id)
                if (any(is.na(j))) stop("filter_data_for_heatmap: Cell mismatch !!!")
                clustCell=clustCell[j,]
                annCellThis$cluster_hclust=clustCell[,paste0("clustId_",nrow(snacsObj$annHash))]
                grp=annCellThis$cluster_hclust
            } else {
                distMat=distfun(t(datThis))
                clustCell=hclust(distMat,method=linkMethod)
                clustCell=cutCluster(clustCell,ann=annCellThis,nClust=nrow(snacsObj$annHash))
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
                        #pvMat[i,gId]=anova(lm(datThis[i,]~grp2))[1,5]
                        res=fisher.test(datThis[i,],grp2)
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
                library(skmeans)
                library(cluster)
                clustCell=skmeans(t(datThis),k=nrow(snacsObj$annHash))
                silCell=silhouette(clustCell)
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
            if (F) {
                j=c()
                for (gId in 1:length(grpUniq)) {
                    jj=which(grp==grpUniq[gId])
                    x=hclust(distfun(t(datThis[,jj])),method=linkMethod)
                    k=x$order
                    j=c(j,jj[k])
                }
            }
            datThis=datThis[,j]
            hashesThis=hashesThis[,j]
            annCellThis=annCellThis[j,]
        }
    }

    ################################################
    if (subsetSnpFlag!="") {
        if (subsetSnpFlag=="_knownHashSnpOrder") {
            tbl=read.table(paste0("../heatmap/_knownHash/clustInfoSnp_knownHash_",snacsObj$exptName,".txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
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
    #datThis[datThis==0]=-1
    #snacsObj=list(mut=datThis,annSNP=annSNPthis,annCell=annCellThis)
    
    names(annCellThis)[match(paste0("cluster_",clustMethod),names(annCellThis))]=paste0("clustRankedSNPs_",clustMethod)

    snacsObj[["mut"]]=datThis
    snacsObj[["hashes"]]=hashesThis
    snacsObj[["annCell"]]=annCellThis
    snacsObj[["annSNP"]]=annSNPthis

    ###########################################################
    ###########################################################
    createHeatmap(snacsObj,cell_anno_var,col_dend=T,row_dend=F,outputfileName=outputfileName)
    
    ###########################################################
    ###########################################################

    invisible(snacsObj)

}


###########################################################
###########################################################
getBestSNPs.internal=function(snacsObj,cell_anno_var,hashColors=NULL,clustMethod=c("hclust","skmean"),pvSnpThres=0) {
    
    #snacsObj=snacsObj; cell_anno_var=c(snacsObj$annHash$hashNames,"clustRankedSNPs_hclust"); pvSnpThres=0.05

    ################################################

    clustMethod=clustMethod[1]
    cellClusterFileName=NA
    subsetSnpFlag=""
    #subsetCellFlag="_allCells"
    subsetCellFlag=""

    #cellClusterFileName=paste0("../heatmap/",subsetCellFlag,"/clustInfoCell",subsetCellFlag,"_pvRanked_",snacsObj$exptName,".txt")
    
    outputfileName=paste0("heatmap",ifelse(clustMethod=="skmean","_skmean",""),"_pvBest_",snacsObj$exptName)
    
    ################################################
    library(skmeans)
    library(RColorBrewer)
    library(marray)
    dirSrc="../src/"
    source(paste0(dirSrc,"heatmap4.R"))
    source(paste0(dirSrc,"heatmapRelated.R"))
    source(paste0(dirSrc,"generate_heatmap.R"))

    library(skmeans)
    getSKmeansDist=function(x) {as.dist(skmeans_xdist(x))}
    distfun=getSKmeansDist
    linkMethod="ward.D2"

    ## Distance functions available in "heatmapRelated.R"
    #getDist=getCosineDist
    #getDist=dist

    ################################################
    datThis=snacsObj$mut
    hashesThis=snacsObj$hashes
    annCellThis=snacsObj$annCell
    annSNPthis=snacsObj$annSNP
    if (F) {
        annCellThis$hash[annCellThis$hash==""]=NA
        annCellThis$hashKnown=annCellThis$hash
        annCellThis$hashKnown[!annCellThis$hash%in%c("CD45-26","CD45-27","CD45-28")]=NA
        annCellThis$hashUnknown=annCellThis$hash
        annCellThis$hashUnknown[annCellThis$hash%in%c("CD45-26","CD45-27","CD45-28")]=NA
    }
    
    ################################################
    if (T) {
        if (is.na(pvSnpThres)) {
            if (clustMethod=="skmean") {
                library(skmeans)
                library(cluster)
                clustCell=skmeans(t(datThis),k=nrow(snacsObj$annHash))
                silCell=silhouette(clustCell)
                clustCell=clustCell$cluster[match(annCellThis$id,names(clustCell$cluster))]
                silCell=silCell[match(annCellThis$id,rownames(silCell)),]
                annCellThis$cluster_skmean=clustCell
                annCellThis$cluster_silhouette=silCell[,"sil_width"]
                grp=annCellThis$cluster_skmean
            } else if (!is.na(cellClusterFileName)) {
                clustCell=read.table(cellClusterFileName, sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=-1)
                j=match(annCellThis$id,clustCell$id)
                if (any(is.na(j))) stop("filter_data_for_heatmap: Cell mismatch !!!")
                clustCell=clustCell[j,]
                annCellThis$cluster_hclust=clustCell[,paste0("clustId_",nrow(snacsObj$annHash))]
                grp=annCellThis$cluster_hclust
            } else {
                distMat=distfun(t(datThis))
                clustCell=hclust(distMat,method=linkMethod)
                clustCell=cutCluster(clustCell,ann=annCellThis,nClust=nrow(snacsObj$annHash))
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
                        #pvMat[i,gId]=anova(lm(datThis[i,]~grp2))[1,5]
                        res=fisher.test(datThis[i,],grp2)
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
                dirResult="../output/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
                dirResult="../output/cellClusterPairPvalueHistogram/"; if (!file.exists(dirResult)) dir.create(file.path(dirResult))
                png(paste0(dirResult,"histogram_cellClusterPair_withRankedSNPs_",snacsObj$exptName,".png"),width=480,height=480)
                x=annSNPthis$minPValue_1VsOtherClusters
                hist(x,freq=T,breaks=length(x)/100,main=paste0(snacsObj$exptName),xlab=paste0("1 cell cluster vs. other ",nrow(snacsObj$annHash)-1,": Minimum p-value"),ylab="Count",cex.main=2,cex.lab=1.5,cex.axis=1.5)
                dev.off()
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
                library(skmeans)
                library(cluster)
                clustCell=skmeans(t(datThis),k=nrow(snacsObj$annHash))
                silCell=silhouette(clustCell)
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
            if (F) {
                j=c()
                for (gId in 1:length(grpUniq)) {
                    jj=which(grp==grpUniq[gId])
                    x=hclust(distfun(t(datThis[,jj])),method=linkMethod)
                    k=x$order
                    j=c(j,jj[k])
                }
            }
            datThis=datThis[,j]
            hashesThis=hashesThis[,j]
            annCellThis=annCellThis[j,]
        }
    }

    ################################################
    if (subsetSnpFlag!="") {
        if (subsetSnpFlag=="_knownHashSnpOrder") {
            tbl=read.table(paste0("../heatmap/_knownHash/clustInfoSnp_knownHash_",snacsObj$exptName,".txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
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
    #snacsObj[["bestSNPs"]]=annSNPthis$id
    i=match(annSNPthis$id,snacsObj$annSNP$id)
    snacsObj$mut=snacsObj$mut[i,]
    snacsObj$annSNP=snacsObj$annSNP[i,]

    ###########################################################
    ###########################################################
    if (F) {
        i=match(snacsObj$bestSNPs,snacsObj$annSNP$id)
        clustObj=createHeatmap(list(mut=snacsObj$mut[i,],annSNP=snacsObj$annSNP[i,],annCell=snacsObj$annCell),cell_anno_var,col_dend=T,row_dend=F,outputfileName)
    }
    clustObj=createHeatmap(snacsObj,cell_anno_var,col_dend=T,row_dend=F,outputfileName=outputfileName)
    snacsObj[["hclustObj_bestSNPs"]]=clustObj$colClust
    

    ###########################################################
    ###########################################################

    invisible(snacsObj)
}


###########################################################
###########################################################
