## Vanessa Kennedy single cell sequencing project

####################################################################
####################################################################
getBestSNPs=function(datObj,col_anno_var,hashColors=NULL,clustMethod="hclust",pvSnpThres=0) {

    # datObj=datObj1; col_anno_var=hashNames; hashColors=c("green3","indianred2","dodgerblue3"); pvSnpThres=pvSnpThres
    ################################################

    sphericalK=NA
    cellClusterFileName=NA
    subsetSnpFlag=""
    subsetCellFlag="_allCells"
    
    #cellClusterFileName=paste0("../heatmap/",subsetCellFlag,"/clustInfoCell",subsetCellFlag,"_pvRanked_",exptName,".txt")
    
    fileNameSuffix=paste0(ifelse(is.na(sphericalK),"","_skmean"),"_pvBest_",exptName)
    
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
    datThis=datObj$mut
    annSNPthis=datObj$annSNP
    annCellThis=datObj$annCell
    annCellThis$hash[annCellThis$hash==""]=NA
    annCellThis$hashKnown=annCellThis$hash
    annCellThis$hashKnown[!annCellThis$hash%in%c("CD45-26","CD45-27","CD45-28")]=NA
    annCellThis$hashUnknown=annCellThis$hash
    annCellThis$hashUnknown[annCellThis$hash%in%c("CD45-26","CD45-27","CD45-28")]=NA
    
    ################################################
    if (T) {
        if (is.na(pvSnpThres)) {
            if (!is.na(sphericalK)) {
                library(skmeans)
                library(cluster)
                clustCell=skmeans(t(datThis),k=3)
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
                annCellThis$cluster_hclust=clustCell$clustId_3
                grp=annCellThis$cluster_hclust
            } else {
                distMat=distfun(t(datThis))
                clustCell=hclust(distMat,method=linkMethod)
                clustCell=cutCluster(clustCell,ann=annCellThis,nClust=3)
                j=match(annCellThis$id,clustCell$id)
                if (any(is.na(j))) stop("filter_data_for_heatmap: Cell mismatch !!!")
                clustCell=clustCell[j,]
                annCellThis$cluster_hclust=clustCell$clustId_3
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
            annSNPthis$maxLOR_2Vs1Cluster=apply(lorMat,1,max,na.rm=T)
            i=rev(order(annSNPthis$minPValue_2Vs1Cluster))
        } else {
            i=which(annSNPthis$minPValue_2Vs1Cluster<=pvSnpThres)
            if (length(i)==0) {
                cat("\n\n No. of SNPs :",nrow(annSNPthis),"\n",sep="")
                cat("SNP p-values:\n")
                print(summary(annSNPthis$minPValue_2Vs1Cluster))
                stop(paste0("No SNPs with p-value <= ",pvSnpThres," !!!"))
            }
        }
        datThis=datThis[i,]
        annSNPthis=annSNPthis[i,]
        
        j=apply(datThis,2,function(x) {mean(x==1)}); j=which(j>0 & j<1) ## Exclude cells with no or all mutations
        datThis=datThis[,j]
        annCellThis=annCellThis[j,]
        
        if (!is.na(sphericalK)) {
            if (!is.na(pvSnpThres)) {
                library(skmeans)
                library(cluster)
                clustCell=skmeans(t(datThis),k=3)
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
            tbl=read.table(paste0("../heatmap/_knownHash/clustInfoSnp_knownHash_",exptName,".txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
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
    datObj=list(mut=datThis,annSNP=annSNPthis,annCell=annCellThis)
    
    save(datObj,file=paste0(dirData,"datObj_for_heatmap",subsetCellFlag,fileNameSuffix,".RData"))

    ###########################################################
    ###########################################################
    cat("best : pvSnpThres ",pvSnpThres,"\n")
    create_heatmap(datObj,col_anno_var,col_dend=T,row_dend=F,hashColors=hashColors,fileNameSuffix,pvSnpThres=pvSnpThres)
    
    ###########################################################
    ###########################################################

    invisible(datObj)
}


###########################################################
###########################################################
