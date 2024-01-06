## SNACS_Manuscript
## Run SNACS

####################################################################
####################################################################
library(SNACS)
library(heatmap4)
library(readxl)
source("../analysis/functions.R")

####################################################################
####################################################################
R.version[['version.string']]
packageVersion("SNACS")

####################################################################
####################################################################

## ------------------------------
## Create SNACSList objects

hashColors <- c("indianred2","green3","dodgerblue3","magenta3")
dirInput <- "../data/data_sequenced/"
fileInput <- dir(dirInput,pattern="csv")
dirOutput <- "../data/"

for (fId in 1:length(fileInput)) {
    if (substr(fileInput[fId],1,6)%in%paste0("SNACS",5:7)) dirDepth <- "../data/Total Depth and Alt Depth with all SNPs/" else dirDepth <- NULL
    snacsObj <- createInitialSNACSobject(fileName=fileInput[fId],dirName=dirInput,hashColors=hashColors,dirDepth=dirDepth)
    print(snacsObj)
    save(snacsObj,file=paste0(dirOutput,"snacsObj_init_",snacsObj$exptName,".RData"))
    rm(snacsObj)
}

## ------------------------------
## Filter data

dirData <- "../data/"
fileList <- paste0("snacsObj_init_SNACS",c(5:7),"_unfilt.RData")
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
fileList <- paste0("snacsObj_filt_SNACS",c(5:7),"_unfilt.RData")
for (fId in 1:length(fileList)) {
    load(paste0(dirData,fileList[fId]))
    cat("\n\n----------------- \n")
    print(snacsObj)
    
    snacsObj <- SNACS::getBestSNPs(snacsObj,backgndThreshold=NA,outputFormat=outputFormat)
    print(snacsObj)
    save(snacsObj,file=paste0(dirData,"snacsObj_getBestSNPs_",snacsObj$exptName,".RData"))
}

## ------------------------------
## Impute SNPs

dirData <- "../data/"
fileList <- paste0("snacsObj_getBestSNPs_SNACS",c(5:7),"_unfilt.RData")
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
fileList <- paste0("snacsObj_imp_SNACS",c(5:7),"_unfilt.RData")
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
fileList <- paste0("snacsObj_clust_SNACS",c(5:7),"_unfilt.RData")
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
fileList <- paste0("snacsObj_SNACS",c(5:7),"_unfilt.RData")
for (fId in 1:length(fileList)) {
    load(paste0(dirData,fileList[fId]))
    cat("\n\n----------------- \n")
    print(snacsObj)
    
    clustObj=SNACS::createHeatmap(snacsObj,cell_anno_var=c("snacsRnd2","snacsRnd1","clustBestSNPs_hclust",snacsObj$annHash$hashNames),cell_anno_name=c("snacsRnd2","snacsRnd1","bestSNPsCluster",snacsObj$annHash$hashNames),col_dend=TRUE,row_dend=FALSE,outputFileName=paste0("heatmap_pvBestWithAnn_",snacsObj$exptName),outputFormat="pdf")
}

## ------------------------------
## Make SNACS + doubletD calls for multi-sample experiments

dirData <- "../data/"
fileList <- paste0("snacsObj_SNACS",c(5:7),"_unfilt.RData")
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
fileList <- paste0("snacsObj_withDoubletD_SNACS",c(5:7),"_unfilt.RData")
for (fId in 1:length(fileList)) {
    load(paste0(dirData,fileList[fId]))
    cat("\n\n----------------- \n")
    print(snacsObj)
    
    snacsObj <- getHTOdemuxCall(snacsObj=snacsObj)
    
    createAnnotatedHeatmap(snacsObj)
}

## ------------------------------
## Get sample color legend

snacsExpt="SNACS7"

dirOutput <- "../output/legend/"

dirData <- "../data/"
fName <- paste0("snacsObj_imp_",snacsExpt,"_unfilt.RData")
load(paste0(dirData,fName))

if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))
pdf(paste0(dirOutput,"sampleColorLegend.pdf"))
heatmap4::sampleColorLegend(tls=snacsObj$annHash$patient,col=snacsObj$annHash$patColors)
dev.off()

####################################################################
####################################################################
