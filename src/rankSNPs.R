## Vanessa Kennedy single cell sequencing project

## Heatmaps: Generated under ../heatmap folder

####################################################################
####################################################################

################################################
rank_SNPs=function(datObj,col_anno_var,hashColors=NULL,clustMethod="hclust") {
    
    pvSnpThres=NA
    sphericalK=NA
    cellClusterFileName=NA
    subsetSnpFlag=""
    subsetCellFlag="_allCells"
    
    fileNameSuffix=paste0(ifelse(is.na(sphericalK),"","_skmean"),"_pvRanked_",exptName)
    
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
    datThis=t(datObj$mutation)
    annThis=datObj$annSNP
    phenThis=datObj$annCell
    phenThis$hash[phenThis$hash==""]=NA
    phenThis$hashKnown=phenThis$hash
    phenThis$hashKnown[!phenThis$hash%in%c("CD45-26","CD45-27","CD45-28")]=NA
    phenThis$hashUnknown=phenThis$hash
    phenThis$hashUnknown[phenThis$hash%in%c("CD45-26","CD45-27","CD45-28")]=NA
    
    ################################################
    if (T) {
        if (is.na(pvSnpThres)) {
            if (!is.na(sphericalK)) {
                library(skmeans)
                library(cluster)
                clustCell=skmeans(t(datThis),k=3)
                silCell=silhouette(clustCell)
                clustCell=clustCell$cluster[match(phenThis$id,names(clustCell$cluster))]
                silCell=silCell[match(phenThis$id,rownames(silCell)),]
                phenThis$cluster_skmean=clustCell
                phenThis$cluster_silhouette=silCell[,"sil_width"]
                grp=phenThis$cluster_skmean
            } else if (!is.na(cellClusterFileName)) {
                clustCell=read.table(cellClusterFileName, sep="\t", h=T, quote="", comment.char="",as.is=T,fill=T,nrow=-1)
                j=match(phenThis$id,clustCell$id)
                if (any(is.na(j))) stop("filter_data_for_heatmap: Cell mismatch !!!")
                clustCell=clustCell[j,]
                phenThis$cluster_hclust=clustCell$clustId_3
                grp=phenThis$cluster_hclust
            } else {
                distMat=distfun(t(datThis))
                clustCell=hclust(distMat,method=linkMethod)
                clustCell=cutCluster(clustCell,ann=phenThis,nClust=3)
                j=match(phenThis$id,clustCell$id)
                if (any(is.na(j))) stop("filter_data_for_heatmap: Cell mismatch !!!")
                clustCell=clustCell[j,]
                phenThis$cluster_hclust=clustCell$clustId_3
                grp=phenThis$cluster_hclust
            }
            grpUniq=unique(grp)
            pvMat=matrix(nrow=nrow(annThis),ncol=length(grpUniq))
            lorMat=matrix(nrow=nrow(annThis),ncol=length(grpUniq))
            for (gId in 1:length(grpUniq)) {
                grp2=rep(1,length(grp))
                grp2[which(!grp%in%grpUniq[gId])]=0
                for (i in 1:nrow(annThis)) {
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
            annThis$minPValue_2Vs1Cluster=apply(pvMat,1,min,na.rm=T)
            annThis$maxLOR_2Vs1Cluster=apply(lorMat,1,max,na.rm=T)
            i=rev(order(annThis$minPValue_2Vs1Cluster))
        } else {
            i=which(annThis$minPValue_2Vs1Cluster<=pvSnpThres)
            if (length(i)==0) {
                cat("\n\n No. of SNPs :",nrow(annThis),"\n",sep="")
                cat("SNP p-values:\n")
                print(summary(annThis$minPValue_2Vs1Cluster))
                stop(paste0("No SNPs with p-value <= ",pvSnpThres," !!!"))
            }
        }
        datThis=datThis[i,]
        annThis=annThis[i,]
        
        j=apply(datThis,2,function(x) {mean(x==1)}); j=which(j>0 & j<1) ## Exclude cells with no or all mutations
        datThis=datThis[,j]
        phenThis=phenThis[j,]
        
        if (!is.na(sphericalK)) {
            if (!is.na(pvSnpThres)) {
                library(skmeans)
                library(cluster)
                clustCell=skmeans(t(datThis),k=3)
                silCell=silhouette(clustCell)
                clustCell=clustCell$cluster[match(phenThis$id,names(clustCell$cluster))]
                silCell=silCell[match(phenThis$id,rownames(silCell)),]
                phenThis$cluster_skmean2=clustCell
                phenThis$cluster_silhouette2=silCell[,"sil_width"]
                grp=phenThis$cluster_skmean2
            }
            grp=phenThis$cluster_skmean
            grpUniq=unique(grp)
            j=c()
            for (gId in 1:length(grpUniq)) {
                jj=which(grp==grpUniq[gId])
                k=order(phenThis$cluster_silhouette[jj])
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
            phenThis=phenThis[j,]
        }
    }

    ################################################
    if (subsetSnpFlag!="") {
        if (subsetSnpFlag=="_knownHashSnpOrder") {
            tbl=read.table(paste0("../heatmap/_knownHash/clustInfoSnp_knownHash_",exptName,".txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            i=match(tbl$id,annThis$id)
            i=rev(i)
            datThis=datThis[i,]
            annThis=annThis[i,]
        }
    }
        
    ################################################
    j=apply(datThis,2,function(x) {mean(x==1)}); j=which(j>0 & j<1) ## Exclude cells with no or all mutations
    datThis=datThis[,j]
    phenThis=phenThis[j,]

    ################################################
    #datThis[datThis==0]=-1
    datObj=list(mut=datThis,ann=annThis,phen=phenThis)
    
    #save(datObj,file=paste0(dirData,"datObj_for_heatmap.RData"))
    save(datObj,file=paste0(dirData,"datObj_for_heatmap",subsetCellFlag,ifelse(is.na(sphericalK),"",paste0("_skmean",sphericalK)),ifelse(is.na(cellClusterFileName),"",paste0("_pvRanked")),ifelse(is.na(pvSnpThres),"",paste0("_pv",pvSnpThres)),fileNameSuffix,"_",exptName,".RData"))
    

    
    ###########################################################
    ###########################################################
    create_heatmap(datObj,col_anno_var,col_dend=T,row_dend=F,hashColors=hashColors,fileNameSuffix)
    
    ###########################################################
    ###########################################################

    invisible(datObj)

}


###########################################################
###########################################################
