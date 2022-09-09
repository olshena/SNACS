####################################################################
####################################################################

createSNACSobject=function(mut,annSNP,annCell,exptName,hashNames,hashColors=NULL) {
    if (is.null(hashColors)) hashColors=rainbow(n=length(hashNames))
    
    if (!is.matrix(mut)) stop("mut has to be a matrix of 0s and 1s")
    if (!any(mut%in%c(0,1))) stop("mut has to be a matrix of 0s and 1s")
    if (nrow(mut)!=nrow(annSNP)) stop("Number of rows in mut and annSNP do not match")
    if (ncol(mut)!=nrow(annCell)) stop("Number of columns in mut and number of rows in annCell do not match")
    if (any(!hashNames%in%names(annCell))) stop("Hash names not found in annCell")
    if (length(hashColors)<length(hashNames)) stop("Number of hash colors should be the same as the number of hashes")
    if ("id"%in%names(annSNP)) warning("Replacing column id in annSNP with generated SNP IDs")
    if ("id"%in%names(annCell)) warning("Replacing column id in annCell with generated cell IDs")
    
    rownames(mut)=annSNP$id=paste0("snp",1:nrow(mut))
    colnames(mut)=annCell$id=paste0("cell",1:ncol(mut))
    if (ncol(annSNP)==1) annSNP$desc=annSNP$id
    if (ncol(annCell)==1) annCell$desc=annCell$id
    for (k in 1:ncol(annSNP)) if (is.factor(annSNP[,k])) annSNP[,k]=as.character(annSNP[,k])
    for (k in 1:ncol(annCell)) if (is.factor(annCell[,k])) annCell[,k]=as.character(annCell[,k])

    hashColors=hashColors[1:length(hashNames)]

    ## Generate scratch folder if it doesn't exist
    dirData="../data/"
    if (!file.exists(dirData)) dir.create(file.path(dirData))

    if (!all(c("id","description")%in%names(annSNP))) {
        stop("annSNP has to be a data frame with at least two coulmns - id and description")
    }

    if (!all(c("id","cell_barcode")%in%names(annCell))) {
        stop("annCell has to be a data frame with at least two coulmns - id and cell_barcode")
    }

    snacsObj=list(mut=mut,annSNP=annSNP,annCell=annCell,exptName=exptName,annHash=data.frame(hashNames=hashNames,hashColors=hashColors,stringsAsFactors=F))
    
    invisible(snacsObj)
}

####################################################################
####################################################################

imputeMissingMutations=function(snacsObj,verbose=F) {
    timeStamp=Sys.time()
    #print(format(timeStamp, "%x %X"))
    
    library(VIM)

    if (verbose) cat("\n\nImputing ",snacsObj$exptName," ...\n",sep="")
    dirData="../data/"
    #load(paste0(dirData,"m.int.",exptName,".RData"))

    ## Generate scratch folder if it doesn't exist
    dirOutput="../output/"
    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))

    numSNP=nrow(snacsObj$mut); numCell=ncol(snacsObj$mut)
    if (verbose) {
        cat("\n Number of SNPs in raw data: ",numSNP,"\n",sep="")
        cat("\n Number of cells in raw data: ",numCell,"\n",sep="")
    }

    propMissPerCell <- apply(snacsObj$mut,2,function(x) sum(is.na(x)))/nrow(snacsObj$mut)
    if (verbose) {
        if (F) {
            cat("\npropMissPerCell")
            print(summary(propMissPerCell))
            #save(propMissPerCell,file=paste0(dirOutput,"propMissPerCell_",snacsObj$exptName,".RData"))
            png(paste0(dirOutput,"hist_propMissPerCell_",snacsObj$exptName,".png"))
            hist(propMissPerCell,xlim=c(0,1),main=snacsObj$exptName,xlab="Proportion missing in a cell",ylab="Count",breaks=100)
            dev.off()
        }
    }
    i=1:nrow(snacsObj$mut); j=which(propMissPerCell<0.4)
    datThis=snacsObj$mut[i,j]; annSNPthis=snacsObj$annSNP[i,]; annCellThis=snacsObj$annCell[j,]
    rm(propMissPerCell)
    
    propMissPerSNP <- apply(datThis,1,function(x) sum(is.na(x)))/ncol(datThis)
    if (verbose) {
        if (F) {
            cat("\npropMissPerSNP")
            print(summary(propMissPerSNP))
            #save(propMissPerSNP,file=paste0(dirOutput,"propMissPerSNP_",snacsObj$exptName,".RData"))
            png(paste0(dirOutput,"propMissPerSNP_",snacsObj$exptName,".png"))
            hist(propMissPerSNP,xlim=c(0,1),main=snacsObj$exptName,xlab="Proportion missing in a SNP",ylab="Count",breaks=100)
            dev.off()
        }
    }
    i=which(propMissPerSNP<0.4); j=1:nrow(annCellThis)
    datThis=datThis[i,j]; annSNPthis=annSNPthis[i,]; annCellThis=annCellThis[j,]
    rm(propMissPerSNP)
    
    #propMutPerSNP <- apply(datThis,1,sum,na.rm=TRUE)/ncol(datThis)
    propMutPerSNP <- apply(datThis,1,function(x) {mean(x==1,na.rm=T)})
    if (verbose) {
        if (F) {
            cat("\npropMutPerSNP")
            print(summary(propMutPerSNP))
        }
    }
    i=which(propMutPerSNP>=0.1 & propMutPerSNP<=0.8); j=1:nrow(annCellThis)
    datThis=datThis[i,j]; annSNPthis=annSNPthis[i,]; annCellThis=annCellThis[j,]
    rm(propMutPerSNP)
    
    if (verbose) {
        cat("\n Number of SNPs after filtering: ",nrow(datThis),"\n",sep="")
        cat("\n Number of cells after filtering: ",ncol(datThis),"\n",sep="")
        #cat("\n Proportion of SNPs kept after filtering: ",round(nrow(datThis)/numSNP,2),"\n",sep="")
        #cat("\n Proportion of cells kept after filtering: ",round(ncol(datThis)/numCell,2),"\n",sep="")
    }
    
    missMat=matrix(F,nrow=nrow(datThis),ncol=ncol(datThis),dimnames=list(rownames(datThis),colnames(datThis))); missMat[is.na(datThis)]=T
    if (T) {
        x=matrix(nrow=nrow(datThis),ncol=ncol(datThis),dimnames=list(rownames(datThis),colnames(datThis)))
        x[datThis==0]="F"; x[datThis==1]="T"
        cellNames=colnames(x)
        rm(datThis)
        x=t(x)
        x=kNN(x,k=5,imp_var=F)
        x=as.matrix(x); x=t(x)
        datThis=matrix(nrow=nrow(x),ncol=ncol(x),dimnames=list(rownames(x),cellNames))
        datThis[x=="F"]=0
        datThis[x=="T"]=1
        rm(x)
    }
    
    #save(datThis,annSNPthis,annCellThis,file=paste0(dirData,"m.int.filt.imp.",snacsObj$exptName,".RData"))
    
    #snacsObj=list(mut=datThis,annSNP=annSNPthis,annCell=annCellThis,missing=missMat)
    snacsObj[["mut"]]=datThis
    snacsObj[["annSNP"]]=annSNPthis
    snacsObj[["annCell"]]=annCellThis
    snacsObj[["missing"]]=missMat

    if (verbose) cat("\nImputation done\n\n",sep="")
    
    timeStamp=c(timeStamp,Sys.time())
    #print(format(timeStamp[2], "%x %X"))
    print(diff(timeStamp))
    
    invisible(snacsObj)
    
}

####################################################################
####################################################################

