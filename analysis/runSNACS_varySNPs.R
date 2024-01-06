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
    ## Filter data

    dirData <- "../data/"
    fileList <- paste0("snacsObj_init_SNACS",c(5:7),exptNameSuffix,".RData")
    for (fId in 1:length(fileList)) {
        load(paste0(dirData,fileList[fId]))
        cat("\n\n----------------- \n")
        print(snacsObj)
        
        snacsObj <- SNACS::filterData(snacsObj=snacsObj)
        print(snacsObj)
        save(snacsObj,file=paste0(dirData,"snacsObj_filt_",snacsObj$exptName,".RData"))
    }
    
    ## ------------------------------
    ## Select best SNPs for multi-sample experiments

    outputFormat <- "pdf"

    dirData <- "../data/"
    fileList <- paste0("snacsObj_filt_SNACS",c(5:7),exptNameSuffix,".RData")
    for (fId in 1:length(fileList)) {
        load(paste0(dirData,fileList[fId]))
        cat("\n\n----------------- \n")
        print(snacsObj)
        
        snacsObj <- SNACS::getBestSNPs(snacsObj,numSNP=numSNP,backgndThreshold=NA,outputFormat=outputFormat)
        print(snacsObj)
        save(snacsObj,file=paste0(dirData,"snacsObj_getBestSNPs_",snacsObj$exptName,".RData"))
    }

    ## ------------------------------
    ## Impute SNPs

    dirData <- "../data/"
    fileList <- paste0("snacsObj_getBestSNPs_SNACS",c(5:7),exptNameSuffix,".RData")
    for (fId in 1:length(fileList)) {
        load(paste0(dirData,fileList[fId]))
        cat("\n\n----------------- \n")
        print(snacsObj)
        
        snacsObj <- SNACS::imputeMissingMutations(snacsObj=snacsObj)
        print(snacsObj)
        save(snacsObj,file=paste0(dirData,"snacsObj_imp_",snacsObj$exptName,".RData"))
    }

    ## ------------------------------
    ## Cluster cells using SNP data

    outputFormat <- "pdf"

    dirData <- "../data/"
    fileList <- paste0("snacsObj_imp_SNACS",c(5:7),exptNameSuffix,".RData")
    for (fId in 1:length(fileList)) {
        load(paste0(dirData,fileList[fId]))
        cat("\n\n----------------- \n")
        print(snacsObj)
        
        snacsObj <- SNACS::clusterCellsWithSNPdata(snacsObj)
        print(snacsObj)
        save(snacsObj,file=paste0(dirData,"snacsObj_clust_",snacsObj$exptName,".RData"))
    }

    ## ------------------------------
    ## Make SNACS calls for multi-sample experiments

    dirData <- "../data/"
    fileList <- paste0("snacsObj_clust_SNACS",c(5:7),exptNameSuffix,".RData")
    for (fId in 1:length(fileList)) {
        load(paste0(dirData,fileList[fId]))
        cat("\n\n----------------- \n")
        print(snacsObj)

        snacsObj <- SNACS::makeSnacsCall(snacsObj)
        print(snacsObj)
        save(snacsObj,file=paste0(dirData,"snacsObj_",snacsObj$exptName,".RData"))
    }

    ## ------------------------------
    ## Generate heatmap with annotation

    outputFormat <- "pdf"
    dirData <- "../data/"
    fileList <- paste0("snacsObj_SNACS",c(5:7),exptNameSuffix,".RData")
    for (fId in 1:length(fileList)) {
        load(paste0(dirData,fileList[fId]))
        cat("\n\n----------------- \n")
        print(snacsObj)
        
        clustObj=SNACS::createHeatmap(snacsObj,cell_anno_var=c("snacsRnd2","snacsRnd1","clustBestSNPs_hclust",snacsObj$annHash$hashNames),cell_anno_name=c("snacsRnd2","snacsRnd1","bestSNPsCluster",snacsObj$annHash$hashNames),col_dend=TRUE,row_dend=FALSE,outputFileName=paste0("heatmap_pvBestWithAnn_",snacsObj$exptName),outputFormat="pdf")
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
