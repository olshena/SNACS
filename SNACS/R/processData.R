####################################################################
####################################################################
#' Creates a SNACSList object.
#'
#' SNACS stores data in a simple list-based data object called a SNACSList.
#'
#' @param mut Matrix having mutation status as 0s and 1s
#' @param hashes Matrix having hash values as numeric data
#' @param exptName Experiment name
#' @param annCell Data frame of cell annotation
#' @param annSNP Data frame of SNP annotation
#' @param annHash Data frame of hash annotation
#' @param hashColors Colors of the hashes
#' @return A SNACSList object
#' @export
SNACSList=function(mut,hashes,exptName,annCell=NULL,annSNP=NULL,annHash=NULL,hashColors=NULL) {
    
    if (!is.matrix(mut)) stop("mut has to be a matrix of 0s and 1s")
    if (!is.numeric(mut[1,1]) | !all(c(0,1)%in%mut) | any(!mut[!is.na(mut)]%in%c(0,1))) stop("mut has to be a matrix of 0s and 1s")
    if (ncol(mut)!=ncol(hashes)) stop("hashes has to be matrix with each row represeting a hash and each column a cells")
    #if ("id"%in%names(annSNP)) warning("Replacing column id in annSNP with generated SNP IDs")
    #if ("id"%in%names(annCell)) warning("Replacing column id in annCell with generated cell IDs")
    
    hashNames=rownames(hashes)
    if (is.null(hashColors)) hashColors=grDevices::rainbow(n=length(hashNames))
    if (length(hashColors)<length(hashNames)) stop("Number of hash colors should be the same as the number of hashes")
    if ("cyan2"%in%hashColors) stop("Cyan is reserved for doublets. Please choose a different hash color")
    if (is.null(annCell)) {
        annCell=data.frame(id=paste0("cell",1:ncol(mut)))
    } else {
        if (ncol(mut)!=nrow(annCell)) stop("Number of columns in mut and number of rows in annCell do not match")
    }
    if (is.null(annSNP)) {
        annSNP=data.frame(id=paste0("snp",1:nrow(mut)))
    } else {
        if (nrow(mut)!=nrow(annSNP)) stop("Number of rows in mut and annSNP do not match")
    }
    colnames(mut)=colnames(hashes)=annCell$id=paste0("cell",1:ncol(mut))
    rownames(mut)=annSNP$id=paste0("snp",1:nrow(mut))
    if (ncol(annCell)==1) annCell$desc=annCell$id
    if (ncol(annSNP)==1) annSNP$desc=annSNP$id
    for (k in 1:ncol(annCell)) if (is.factor(annCell[,k])) annCell[,k]=as.character(annCell[,k])
    for (k in 1:ncol(annSNP)) if (is.factor(annSNP[,k])) annSNP[,k]=as.character(annSNP[,k])

    hashColors=hashColors[1:length(hashNames)]
    tbl=data.frame(hashNames=hashNames,hashColors=hashColors,stringsAsFactors=F)
    if (is.null(annHash)) {
        annHash=tbl
    } else {
        if (length(hashNames)!=nrow(annHash)) stop("Number of hashes in hashes matrix and annHash do not match")
        k=which(!names(annHash)%in%names(tbl))
        if (length(k)==0) {
            annHash=tbl
        } else {
            nm=c(names(tbl),names(annHash)[k])
            annHash=cbind(tbl,annHash[,k],stringsAsFactors=F)
            names(annHash)=nm
            rm(nm)
        }
    }
    rm(tbl)
    
    ## Generate scratch folder if it doesn't exist
    dirData="../data/"
    if (!file.exists(dirData)) dir.create(file.path(dirData))

    if (!all(c("id","desc")%in%names(annSNP))) {
        stop("annSNP has to be a data frame with at least two columns - id and description")
    }

    if (!all(c("id","desc")%in%names(annCell))) {
        stop("annCell has to be a data frame with at least two columns - id and desc")
    }

    snacsObj=list(mut=mut,hashes=hashes,exptName=exptName,annHash=annHash,annCell=annCell,annSNP=annSNP)
    
    #snacsObj=new("SNACSList",snacsObj)
    attr(snacsObj,"class")="SNACSList"

    invisible(snacsObj)
}

####################################################################
####################################################################
#' Prints SNACSList object.
#'
#' @method print SNACSList
#' @param x SNACSList object
#' @param ... ignored
#' @export
print.SNACSList <- function(x,...) {
    cat('An object of class "SNACSList"\n',sep="")
    cat("Experiment name: ",x$exptName, "\n",sep="")
    cat("No. of SNPs: ",nrow(x$mut),"\n",sep="")
    cat("No. of cells: ",ncol(x$mut),"\n",sep="")
    cat("Hashes: ",paste(x$annHash$hashNames,collapse=", "),"\n",sep="")
    cat("SNACSList attributes: ",paste(names(x),collapse=", "),"\n",sep="")
    cat('Output files are saved in "../output" folder\n',sep="")
}

####################################################################
####################################################################
#' Imputes missing mutations in SNACSList object.
#'
#' @param snacsObj SNACSList object
#' @param verbose Prints information when running the method
#' @param proportionMissingPerCell .
#' @param proportionMissingPerSNP .
#' @param proportionMutatedPerCell .
#' @param proportionMutatedPerSNP .
#' @return A SNACSList object
#' @export
imputeMissingMutations=function(snacsObj,proportionMissingPerCell=0.4,proportionMissingPerSNP=0.4,proportionMutatedPerCell=c(0,1),proportionMutatedPerSNP=c(0.1,0.8),verbose=F) {
    timeStamp=Sys.time()
    #print(format(timeStamp, "%x %X"))

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

    propMissPerCellVec <- apply(snacsObj$mut,2,function(x) sum(is.na(x)))/nrow(snacsObj$mut)
    if (verbose) {
        if (F) {
            cat("\npropMissPerCellVec")
            print(summary(propMissPerCellVec))
            #save(propMissPerCellVec,file=paste0(dirOutput,"propMissPerCellVec_",snacsObj$exptName,".RData"))
            grDevices::png(paste0(dirOutput,"hist_propMissPerCellVec_",snacsObj$exptName,".png"))
            graphics::hist(propMissPerCellVec,xlim=c(0,1),main=snacsObj$exptName,xlab="Proportion missing in a cell",ylab="Count",breaks=100)
            grDevices::dev.off()
        }
    }
    i=1:nrow(snacsObj$mut); j=which(propMissPerCellVec<proportionMissingPerCell)
    if (length(j)<10) {
        cat("\nProportion missing per cell:\n")
        print(round(summary(propMissPerCellVec),2))
        stop(paste0("Not enough cells having proportion missing below ",proportionMissingPerCell))
    }
    datThis=snacsObj$mut[i,j]; hashesThis=snacsObj$hashes[,j]; annCellThis=snacsObj$annCell[j,]; annSNPthis=snacsObj$annSNP[i,]
    rm(propMissPerCellVec)
    
    propMissPerSNPvec <- apply(datThis,1,function(x) sum(is.na(x)))/ncol(datThis)
    if (verbose) {
        if (F) {
            cat("\npropMissPerSNPvec")
            print(summary(propMissPerSNPvec))
            #save(propMissPerSNPvec,file=paste0(dirOutput,"propMissPerSNPvec_",snacsObj$exptName,".RData"))
            grDevices::png(paste0(dirOutput,"propMissPerSNPvec_",snacsObj$exptName,".png"))
            graphics::hist(propMissPerSNPvec,xlim=c(0,1),main=snacsObj$exptName,xlab="Proportion missing in a SNP",ylab="Count",breaks=100)
            grDevices::dev.off()
        }
    }
    i=which(propMissPerSNPvec<proportionMissingPerSNP); j=1:nrow(annCellThis)
    if (length(i)<10) {
        cat("\nProportion missing per SNP:\n")
        print(round(summary(proportionMissingPerSNP),2))
        stop(paste0("Not enough SNPs having proportion missing below ",proportionMissingPerSNP))
    }
    datThis=datThis[i,j]; hashesThis=hashesThis[,j]; annCellThis=annCellThis[j,]; annSNPthis=annSNPthis[i,]
    rm(propMissPerSNPvec)
    
    #propMutPerSNPvec <- apply(datThis,1,sum,na.rm=TRUE)/ncol(datThis)
    propMutPerSNPvec <- apply(datThis,1,function(x) {mean(x==1,na.rm=T)})
    if (verbose) {
        if (F) {
            cat("\npropMutPerSNPvec")
            print(summary(propMutPerSNPvec))
        }
    }
    i=which(propMutPerSNPvec>=proportionMutatedPerSNP[1] & propMutPerSNPvec<=proportionMutatedPerSNP[2]); j=1:nrow(annCellThis)
    if (length(i)<10) {
        cat("\nProportion mutated per SNP:\n")
        print(round(summary(propMutPerSNPvec),2))
        stop(paste0("Not enough SNPs having proportion mutated in the range ",paste(proportionMutatedPerSNP,collapse=" - ")))
    }
    datThis=datThis[i,j]; hashesThis=hashesThis[,j]; annCellThis=annCellThis[j,]; annSNPthis=annSNPthis[i,]
    rm(propMutPerSNPvec)
    
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
        x=VIM::kNN(x,k=5,imp_var=F)
        x=as.matrix(x); x=t(x)
        datThis=matrix(nrow=nrow(x),ncol=ncol(x),dimnames=list(rownames(x),cellNames))
        datThis[x=="F"]=0
        datThis[x=="T"]=1
        rm(x)
    }

    propMutPerCellVec=apply(datThis,2,function(x) mean(x==1,na.rm=T)); j=which(propMutPerCellVec>proportionMutatedPerCell[1] & propMutPerCellVec<proportionMutatedPerCell[2])
    if (length(i)<10) {
        cat("\nProportion mutated per cell:\n")
        print(round(summary(propMutPerCellVec),2))
        stop(paste0("Not enough cells having proportion mutated in the range ",paste(proportionMutatedPerCell,collapse=" - ")))
    }
    datThis=datThis[,j]; hashesThis=hashesThis[,j]; annCellThis=annCellThis[j,]
    rm(propMutPerCellVec)

    snacsObj[["mut"]]=datThis
    snacsObj[["hashes"]]=hashesThis
    snacsObj[["annCell"]]=annCellThis
    snacsObj[["annSNP"]]=annSNPthis
    snacsObj[["missing"]]=missMat

    if (verbose) {
        cat("\nImputation done\n\n",sep="")
        
        timeStamp=c(timeStamp,Sys.time())
        #print(format(timeStamp[2], "%x %X"))
        print(diff(timeStamp))
    }
    
    invisible(snacsObj)
    
}

####################################################################
####################################################################

