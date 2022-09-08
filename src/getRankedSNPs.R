getRankedSNPs=function(snacsObj,col_anno_var,clustMethod=c("hclust","skmean")) {
    
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
    annCellThis=annCellThis[j,]

    ################################################
    #datThis[datThis==0]=-1
    #snacsObj=list(mut=datThis,annSNP=annSNPthis,annCell=annCellThis)
    
    names(annCellThis)[match(paste0("cluster_",clustMethod),names(annCellThis))]=paste0("clustRankedSNPs_",clustMethod)

    snacsObj[["mut"]]=datThis
    snacsObj[["annSNP"]]=annSNPthis
    snacsObj[["annCell"]]=annCellThis
    
    ###########################################################
    ###########################################################
    createHeatmap(snacsObj,col_anno_var,col_dend=T,row_dend=F,outputfileName=outputfileName)
    
    ###########################################################
    ###########################################################

    invisible(snacsObj)

}


###########################################################
###########################################################
