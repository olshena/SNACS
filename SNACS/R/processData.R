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
#' @param hashColors Colors of the hashes
#' @return A SNACSList object
#' @export
SNACSList=function(mut,hashes,exptName,annCell=NULL,annSNP=NULL,hashColors=NULL) {
    
    if (!is.matrix(mut)) stop("mut has to be a matrix of 0s and 1s")
    if (!any(mut%in%c(0,1))) stop("mut has to be a matrix of 0s and 1s")
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

    ## Generate scratch folder if it doesn't exist
    dirData="../data/"
    if (!file.exists(dirData)) dir.create(file.path(dirData))

    if (!all(c("id","desc")%in%names(annSNP))) {
        stop("annSNP has to be a data frame with at least two columns - id and description")
    }

    if (!all(c("id","desc")%in%names(annCell))) {
        stop("annCell has to be a data frame with at least two columns - id and desc")
    }

    snacsObj=list(mut=mut,hashes=hashes,exptName=exptName,annHash=data.frame(hashNames=hashNames,hashColors=hashColors,stringsAsFactors=F),annCell=annCell,annSNP=annSNP)
    
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
#' @return A SNACSList object
#' @export
imputeMissingMutations=function(snacsObj,verbose=F) {
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

    propMissPerCell <- apply(snacsObj$mut,2,function(x) sum(is.na(x)))/nrow(snacsObj$mut)
    if (verbose) {
        if (F) {
            cat("\npropMissPerCell")
            print(summary(propMissPerCell))
            #save(propMissPerCell,file=paste0(dirOutput,"propMissPerCell_",snacsObj$exptName,".RData"))
            grDevices::png(paste0(dirOutput,"hist_propMissPerCell_",snacsObj$exptName,".png"))
            graphics::hist(propMissPerCell,xlim=c(0,1),main=snacsObj$exptName,xlab="Proportion missing in a cell",ylab="Count",breaks=100)
            grDevices::dev.off()
        }
    }
    i=1:nrow(snacsObj$mut); j=which(propMissPerCell<0.4)
    datThis=snacsObj$mut[i,j]; hashesThis=snacsObj$hashes[,j]; annCellThis=snacsObj$annCell[j,]; annSNPthis=snacsObj$annSNP[i,]
    rm(propMissPerCell)
    
    propMissPerSNP <- apply(datThis,1,function(x) sum(is.na(x)))/ncol(datThis)
    if (verbose) {
        if (F) {
            cat("\npropMissPerSNP")
            print(summary(propMissPerSNP))
            #save(propMissPerSNP,file=paste0(dirOutput,"propMissPerSNP_",snacsObj$exptName,".RData"))
            grDevices::png(paste0(dirOutput,"propMissPerSNP_",snacsObj$exptName,".png"))
            graphics::hist(propMissPerSNP,xlim=c(0,1),main=snacsObj$exptName,xlab="Proportion missing in a SNP",ylab="Count",breaks=100)
            grDevices::dev.off()
        }
    }
    i=which(propMissPerSNP<0.4); j=1:nrow(annCellThis)
    datThis=datThis[i,j]; hashesThis=hashesThis[,j]; annCellThis=annCellThis[j,]; annSNPthis=annSNPthis[i,]
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
    datThis=datThis[i,j]; hashesThis=hashesThis[,j]; annCellThis=annCellThis[j,]; annSNPthis=annSNPthis[i,]
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
        x=VIM::kNN(x,k=5,imp_var=F)
        x=as.matrix(x); x=t(x)
        datThis=matrix(nrow=nrow(x),ncol=ncol(x),dimnames=list(rownames(x),cellNames))
        datThis[x=="F"]=0
        datThis[x=="T"]=1
        rm(x)
    }

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

