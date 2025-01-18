## SNACS_Manuscript
## Run SNACS on the core multi-sample experiments (Experiments 5-7) described in section 3.1 of the manuscript

####################################################################
####################################################################
library(SNACS)
#library(heatmap4)
#library(readxl)

#source("../../SNACS/analysis/functions.R")
#source("../analysis/functions.R")

####################################################################
####################################################################
R.version[['version.string']]
packageVersion("SNACS")

####################################################################
####################################################################

## ------------------------------
## Extract raw data from hdf5 file. Then create SNACSList objects

hashColors <- c("indianred2","green3","dodgerblue3","magenta3")
dirDataRaw <- "../data/data_sequenced/"
dirData <- "../data/"

for (exptName in paste0("SNACS",5:7)) {
    cat("\n\n------------- ",exptName,"\n")
    switch(exptName,
        "SNACS1"={
            fileName="GSM8066757.hdf5"
            hashNames=c("TS.1","TS.2")
        },
        "SNACS2"={
            fileName="GSM8066758.hdf5"
            hashNames=c("TS.2","TS.3")
        },
        "SNACS3"={
            fileName="GSM8066759.hdf5"
            hashNames=c("TS.3","TS.4")
        },
        "SNACS4"={
            fileName="GSM8066760.hdf5"
            hashNames=c("TS.1","TS.4")
        },
        "SNACS5"={
            fileName="GSM8066761.hdf5"
            hashNames=c("TS.1","TS.2")
        },
        "SNACS6"={
            fileName="GSM8066762.hdf5"
            hashNames=c("TS.2","TS.3","TS.4")
        },
        "SNACS7"={
            fileName="GSM8066763.hdf5"
            hashNames=c("TS.1","TS.2","TS.3","TS.4")
        }
    )
    h5toList=h5readForSNACS(file=paste0(dirDataRaw,fileName))
    
    snacsObj=SNACSList(mut=h5toList$mut,hashes=h5toList$hashes[hashNames,],exptName=exptName,hashColors=hashColors[match(hashNames,c("TS.1","TS.2","TS.3","TS.4"))],
                          depthTotal=h5toList$depthTotal,depthAlt=h5toList$depthAlt,annCell=h5toList$annCell,annSNP=h5toList$annSNP)
    rm(h5toList)
    snacsObj$annSNP$desc=getSNPdesc(snacsObj$annSNP$desc)
    save(snacsObj,file=paste0(dirData,"snacsObj_init_",exptName,"_unfilt.RData"))
}

## ------------------------------
## Run SNACS

dirData <- "../data/"
fileList <- paste0("snacsObj_init_SNACS",c(5:7),"_unfilt.RData")
for (fId in 1:length(fileList)) {
    load(paste0(dirData,fileList[fId]))
    cat("\n\n----------------- \n")
    print(snacsObj)
    
    snacsObj <- SNACS::runSNACS(snacsObj=snacsObj)
    print(snacsObj)
    save(snacsObj,file=paste0(dirData,"snacsObj_",snacsObj$exptName,".RData"))
}

## ------------------------------
## Generate heatmap with cell annotation

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
## Generate heatmap with SNACS + doubletD annotation

outputFormat <- "pdf"
dirData <- "../data/"
fileList <- paste0("snacsObj_withDoubletD_SNACS",c(5:7),"_unfilt.RData")
for (fId in 1:length(fileList)) {
    load(paste0(dirData,fileList[fId]))
    cat("\n\n----------------- \n")
    print(snacsObj)
    
    clustObj=SNACS::createHeatmap(snacsObj,cell_anno_var=c("snacsPlusDoubletD","doubletD","snacsRnd2","snacsRnd1","clustBestSNPs_hclust",snacsObj$annHash$hashNames),cell_anno_name=c("snacs+DD","doubletD","snacsRnd2","snacsRnd1","bestSNPsCluster",snacsObj$annHash$hashNames),col_dend=TRUE,row_dend=FALSE,outputFileName=paste0("heatmap_",snacsObj$exptName),outputFormat="pdf")
}

## ------------------------------
## Get sample color legend

snacsExpt="SNACS7"

dirOutput <- "../output/legend/"

dirData <- "../data/"
fName <- paste0("snacsObj_",snacsExpt,"_unfilt.RData")
load(paste0(dirData,fName))

if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))
pdf(paste0(dirOutput,"sampleColorLegend.pdf"))
heatmap4::sampleColorLegend(tls=snacsObj$annHash$patient,col=snacsObj$annHash$patColors)
dev.off()

####################################################################
####################################################################
