####################################################################
####################################################################
#' Creates a SNACSList object.
#'
#' SNACS stores data in a simple list-based data object called a SNACSList.
#'
#' @param mut Integer matrix. Matrix having mutation status as 0s and 1s. Rows depict SNPs and columns  cells
#' @param hashes Numeric matrix. Matrix having hash values as numeric data. Rows depict hashes and columns cells
#' @param exptName Character. Experiment name
#' @param annCell Data frame. Table of cell annotation. Default is NULL
#' @param annSNP Data frame. Table of SNP annotation. Default is NULL
#' @param annHash Data frame. Table of hash annotation. Default is NULL
#' @param hashColors Colors of the hashes. Default is NULL
#' @return A SNACSList object
#' @export
SNACSList=function(mut,hashes,exptName,annCell=NULL,annSNP=NULL,annHash=NULL,hashColors=NULL) {
    
    if (!is.matrix(mut)) stop("mut has to be a matrix of 0s and 1s")
    if (!is.numeric(mut[1,1]) | !all(c(0,1)%in%mut) | any(!mut[!is.na(mut)]%in%c(0,1))) stop("mut has to be a matrix of 0s and 1s")
    if (ncol(mut)!=ncol(hashes)) stop("hashes has to be matrix with each row represeting a hash and each column a cells")
    if ("id"%in%names(annSNP)) warning("Replacing column id in annSNP with generated SNP IDs")
    if ("id"%in%names(annCell)) warning("Replacing column id in annCell with generated cell IDs")
    
    hashNames=rownames(hashes)
    if (is.null(hashColors)) hashColors=grDevices::rainbow(n=length(hashNames))
    if (length(hashColors)<length(hashNames)) stop("Number of hash colors should be the same as the number of hashes")
    #if ("cyan2"%in%hashColors) stop("Cyan is reserved for multiplets. Please choose a different hash color")
    cat("Shades of gray are reserved for multiplets\n")
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
    if (!is.character(annCell$id)) annCell$id=as.character(annCell$id)
    if (!is.character(annSNP$id)) annSNP$id=as.character(annSNP$id)

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
    
    ## Generate data folder if it doesn't exist
    dirData="../data/"
    if (!file.exists(dirData)) dir.create(file.path(dirData))

    ## Generate output folder if it doesn't exist
    dirOutput="../output/"
    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))

    if (!all(c("id","desc")%in%names(annSNP))) {
        stop("annSNP has to be a data frame with at least two columns - id and desc")
    }

    if (!all(c("id","desc")%in%names(annCell))) {
        stop("annCell has to be a data frame with at least two columns - id and desc")
    }

    snacsObj=list(mut=mut,hashes=hashes,exptName=exptName,annHash=annHash,annCell=annCell,annSNP=annSNP,processLevels="raw")
    
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
print.SNACSList=function(x,...) {
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
#' Filters mutation data in SNACSList object.
#'
#' @param snacsObj SNACSList object
#' @param proportionMissingPerCell Numeric. Only cells with lower than this proportion of missing values are considered. Default is 0.4. Range is 0-1
#' @param proportionMissingPerSNP Numeric. Only SNPs with lower than this proportion of missing values are considered. Default is 0.4. Range is 0-1
#' @param proportionMutatedPerCell Numeric vector. Only cells having proportion of mutations within this range are considered. Default is c(0,1). Range is 0-1
#' @param proportionMutatedPerSNP Numeric vector. Only SNPs having proportion of mutations within this range are considered. Default is c(0.1,0.8). Range is 0-1
#' @param verbose Logical. Prints information when running the method. Default is FALSE
#' @return A SNACSList object
#' @export
filterData=function(snacsObj,proportionMissingPerCell=0.4,proportionMissingPerSNP=0.4,proportionMutatedPerCell=c(0,1),proportionMutatedPerSNP=c(0.1,0.8),verbose=FALSE) {
    if (verbose) cat("\n\nImputing ",snacsObj$exptName," ...\n",sep="")
    dirData="../data/"

    numSNP=nrow(snacsObj$mut); numCell=ncol(snacsObj$mut)
    if (verbose) {
        cat("\n Number of SNPs in raw data: ",numSNP,"\n",sep="")
        cat("\n Number of cells in raw data: ",numCell,"\n",sep="")
    }

    propMissPerCellVec=apply(snacsObj$mut,2,function(x) sum(is.na(x)))/nrow(snacsObj$mut)
    i=1:nrow(snacsObj$mut); j=which(propMissPerCellVec<proportionMissingPerCell)
    if (length(j)<10) {
        cat("\nProportion missing per cell:\n")
        print(round(summary(propMissPerCellVec),2))
        stop(paste0("Not enough cells having proportion missing below ",proportionMissingPerCell))
    }
    datThis=snacsObj$mut[i,j]; hashesThis=snacsObj$hashes[,j]; annCellThis=snacsObj$annCell[j,]; annSNPthis=snacsObj$annSNP[i,]
    rm(propMissPerCellVec)
    
    propMissPerSNPvec=apply(datThis,1,function(x) sum(is.na(x)))/ncol(datThis)
    i=which(propMissPerSNPvec<proportionMissingPerSNP); j=1:nrow(annCellThis)
    if (length(i)<10) {
        cat("\nProportion missing per SNP:\n")
        print(round(summary(proportionMissingPerSNP),2))
        stop(paste0("Not enough SNPs having proportion missing below ",proportionMissingPerSNP))
    }
    datThis=datThis[i,j]; hashesThis=hashesThis[,j]; annCellThis=annCellThis[j,]; annSNPthis=annSNPthis[i,]
    rm(propMissPerSNPvec)
    
    propMutPerSNPvec=apply(datThis,1,function(x) {mean(x==1,na.rm=T)})
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
    }

    propMutPerCellVec=apply(datThis,2,function(x) mean(x==1,na.rm=T)); j=which(propMutPerCellVec>proportionMutatedPerCell[1] & propMutPerCellVec<proportionMutatedPerCell[2])
    if (length(j)<10) {
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
    snacsObj[["processLevels"]]=c(snacsObj[["processLevels"]],"filtered")
    
    invisible(snacsObj)
}

####################################################################
####################################################################
#' Imputes missing mutations in SNACSList object.
#'
#' @param snacsObj SNACSList object
#' @param verbose Prints information when running the method
#' @return A SNACSList object
#' @export
imputeMissingMutations=function(snacsObj,verbose=F) {
    timeStamp=Sys.time()

    if (!is.na(match("filtered",snacsObj[["processLevels"]]))) stop("Run filterData() before imputing data\n")

    if (verbose) cat("\n\nImputing ",snacsObj$exptName," ...\n",sep="")
    dirData="../data/"
    
    datThis=snacsObj$mut; hashesThis=snacsObj$hashes; annCellThis=snacsObj$annCell; annSNPthis=snacsObj$annSNP

    missMat=matrix(F,nrow=nrow(datThis),ncol=ncol(datThis),dimnames=list(rownames(datThis),colnames(datThis))); missMat[is.na(datThis)]=T
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

    snacsObj[["mut"]]=datThis
    snacsObj[["hashes"]]=hashesThis
    snacsObj[["annCell"]]=annCellThis
    snacsObj[["annSNP"]]=annSNPthis
    snacsObj[["missing"]]=missMat

    if (verbose) {
        cat("\nImputation done\n\n",sep="")
        timeStamp=c(timeStamp,Sys.time())
        print(diff(timeStamp))
    }
    
    invisible(snacsObj)
    
}

####################################################################
####################################################################

