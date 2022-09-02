## Vanessa Kennedy single cell sequencing project

## Process and analyze data
## mpalDat_genotypes.csv should be in ../data folder

####################################################################
####################################################################

process_data=function(mut,annSNP,annCell) {

#inputFile=paste0(mpalName,"_genotypes.csv"); cellAnnFile=paste0(mpalName,"_genotypes_raw.csv"); cellAnnVar=c("CD45.26","CD45.27","CD45.28"); snpAnnFile=paste0("variant_info_",mpalName,".csv"); depthFile=paste0("depth_",mpalName,".csv"); qualityFile=paste0("quality_",mpalName,".csv"); outputFileName=mpalName
#mpalName="giltven1_run3"; inputFile=paste0(mpalName,".csv"); cellAnnFile=""; cellAnnVar=c("TS.1","TS.2"); snpAnnFile=""; depthFile=""; qualityFile=""; outputFileName=mpalName


    ## Generate scratch folder if it doesn't exist
    dirData="../data/"
    if (!file.exists(dirData)) dir.create(file.path(dirData))


    if (!all(c("id","description")%in%names(annSNP))) {
        stop("annSNP has to be a data frame with at least two coulmns - id and description")
    }

    if (!all(c("id","cell_barcode")%in%names(annCell))) {
        stop("annCell has to be a data frame with at least two coulmns - id and cell_barcode")
    }
    
    colnames(mut)=annSNP$id
    rownames(mut)=annCell$id
    
    mut=t(mut)

    datObj=list(mut=mut,annSNP=annSNP,annCell=annCell)
    
    invisible(datObj)
}

####################################################################
####################################################################

impute_data=function(outputFileName="mpal3") {
    library(VIM)
    
    cat("\n\n---------- ",outputFileName," ----------\n",sep="")
    dirData="../data/"; outputFileName=mpalName
    load(paste0(dirData,"m.int.",outputFileName,".RData"))

    ## Generate scratch folder if it doesn't exist
    dirOutput="../output/"
    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))

    numSNP=ncol(datObj$mut); numCell=nrow(datObj$mut)
    
    propMissPerCell <- apply(datObj$mut,1,function(x) sum(is.na(x)))/ncol(datObj$mut)
    cat("\npropMissPerCell")
    print(summary(propMissPerCell))
    save(propMissPerCell,file=paste0(dirOutput,"propMissPerCell_",outputFileName,".RData"))
    png(paste0(dirOutput,"hist_propMissPerCell_",outputFileName,".png"))
    hist(propMissPerCell,xlim=c(0,1),main=outputFileName,xlab="Proportion missing in a cell",ylab="Count",breaks=100)
    dev.off()
    i=1:nrow(annSNP); j=which(propMissPerCell<0.4)
    datObj$mut=datObj$mut[j,i]; annSNP=datObj$annSNP[i,]; annCell=datObj$annCell[j,]
    rm(propMissPerCell)
    
    propMissPerSNP <- apply(datObj$mut,2,function(x) sum(is.na(x)))/nrow(datObj$mut)
    cat("\npropMissPerSNP")
    print(summary(propMissPerSNP))
    save(propMissPerSNP,file=paste0(dirOutput,"propMissPerSNP_",outputFileName,".RData"))
    png(paste0(dirOutput,"propMissPerSNP_",outputFileName,".png"))
    hist(propMissPerSNP,xlim=c(0,1),main=outputFileName,xlab="Proportion missing in a SNP",ylab="Count",breaks=100)
    dev.off()
    i=which(propMissPerSNP<0.4); j=1:nrow(annCell)
    mmpalDat.int=datObj$mut[j,i]; annSNP=annSNP[i,]; annCell=annCell[j,]
    rm(propMissPerSNP)
    
    #propMutPerSNP <- apply(mmpalDat.int,2,sum,na.rm=TRUE)/nrow(mmpalDat.int)
    propMutPerSNP <- apply(mmpalDat.int,2,function(x) {mean(x==1,na.rm=T)})
    cat("\npropMutPerSNP")
    print(summary(propMutPerSNP))
    i=which(propMutPerSNP>=0.1 & propMutPerSNP<=0.8); j=1:nrow(annCell)
    mmpalDat.int=mmpalDat.int[j,i]; annSNP=annSNP[i,]; annCell=annCell[j,]
    rm(propMutPerSNP)
    
    cat("\n no. of SNPs ",round(ncol(mmpalDat.int),2),"\n",sep="")
    cat("\n no. of cells ",round(nrow(mmpalDat.int),2),"\n",sep="")
    cat("\n prop SNPs ",round(ncol(mmpalDat.int)/numSNP,2),"\n",sep="")
    cat("\n prop cells ",round(nrow(mmpalDat.int)/numCell,2),"\n",sep="")
    
    if (F) {
        if (imputeFlag=="snp") {
            mmpalDat.int.filt.imp=as.data.frame(mmpalDat.int)
            rm(mmpalDat.int)
        } else {
            mmpalDat.int.filt.imp=t(mmpalDat.int)
            rm(mmpalDat.int)
            mmpalDat.int.filt.imp=as.data.frame(mmpalDat.int.filt.imp)
        }
        for (k in 1:ncol(mmpalDat.int.filt.imp)) {
            mmpalDat.int.filt.imp[,k]=as.factor(mmpalDat.int.filt.imp[,k])
        }
        mmpalDat.int.filt.imp=knnImputation2(mmpalDat.int.filt.imp,k=10)
        for (k in 1:ncol(mmpalDat.int.filt.imp)) {
            mmpalDat.int.filt.imp[,k]=as.integer(as.character(mmpalDat.int.filt.imp[,k]))
        }
        if (imputeFlag=="snp") {
            mmpalDat.int.filt.imp=as.matrix(mmpalDat.int.filt.imp)
        } else {
            mmpalDat.int.filt.imp=as.matrix(mmpalDat.int.filt.imp)
            mmpalDat.int.filt.imp=t(mmpalDat.int.filt.imp)
        }
    }
    
    x=matrix(nrow=nrow(mmpalDat.int),ncol=ncol(mmpalDat.int),dimnames=list(rownames(mmpalDat.int),colnames(mmpalDat.int)))
    x[mmpalDat.int==0]="F"; x[mmpalDat.int==1]="T"
    rowNames=rownames(x)
    rm(mmpalDat.int)
    x=kNN(x,k=5,imp_var=F)
    x=as.matrix(x)
    mmpalDat.int.filt.imp=matrix(nrow=nrow(x),ncol=ncol(x),dimnames=list(rowNames,colnames(x)))
    mmpalDat.int.filt.imp[x=="F"]=0
    mmpalDat.int.filt.imp[x=="T"]=1
    
    save(mmpalDat.int.filt.imp,annSNP,annCell,file=paste0(dirData,"m.int.filt.imp.",outputFileName,".RData"))
    
    datObj=list(mut=mmpalDat.int.filt.imp,annSNP=annSNP,annCell=annCell)
    
    invisible(datObj)
    
}

####################################################################
####################################################################

