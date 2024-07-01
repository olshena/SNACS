## SNACS_Manuscript
## Vary number of SNPs

####################################################################
####################################################################
library(SNACS)
library(heatmap4)
library(readxl)
source("../analysis/functions.R")

####################################################################
####################################################################
## Parameter

numSNPvec=1:5

for (numSNP in numSNPvec) {
    exptNameSuffix=paste0("_hashFiltNumSNP",numSNP)
    ## ------------------------------
    ## Run SNACS

    dirData <- "../data/"
    fileList <- paste0("snacsObj_init_SNACS",c(5:7),"_unfilt.RData")
    for (eId in c(5:7)) {
        load(paste0(dirData,paste0("snacsObj_init_SNACS",eId,"_unfilt.RData")))
        snacsObj$exptName=paste0("SNACS",eId,exptNameSuffix)
        save(snacsObj,file=paste0("snacsObj_init_SNACS",eId,exptNameSuffix,".RData")
        
        cat("\n\n----------------- \n")
        print(snacsObj)
        
        snacsObj <- SNACS::runSNACS(snacsObj=snacsObj)
        print(snacsObj)
        save(snacsObj,file=paste0(dirData,"snacsObj_",snacsObj$exptName,".RData"))
    }

    ## ------------------------------
    ## Make SNACS + doubletD calls for multi-sample experiments

    dirData <- "../data/"
    fileList <- paste0("snacsObj_SNACS",c(5:7),exptNameSuffix,".RData")
    for (fId in 1:length(fileList)) {
        load(paste0(dirData,fileList[fId]))
        cat("\n\n----------------- \n")
        print(snacsObj)

        snacsObj <- SNACS::runSNACSplusDoubletD(snacsObj)
        print(snacsObj)
        save(snacsObj,file=paste0(dirData,"snacsObj_withDoubletD_",snacsObj$exptName,".RData"))
    }

    ## ------------------------------
    ## Generate heatmap with annotation

    dirData <- "../data/"
    fileList <- paste0("snacsObj_withDoubletD_SNACS",c(5:7),exptNameSuffix,".RData")
    for (fId in 1:length(fileList)) {
        load(paste0(dirData,fileList[fId]))
        cat("\n\n----------------- \n")
        print(snacsObj)
        
        snacsObj <- getHTOdemuxCall(snacsObj=snacsObj)
        
        createAnnotatedHeatmap(snacsObj)
    }
}

####################################################################
####################################################################
