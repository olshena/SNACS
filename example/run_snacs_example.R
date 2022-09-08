
####################################################################
####################################################################
source("../src/heatmapRelated.R")
source("../src/processData.R")
source("../src/getRankedSNPs.R")
source("../src/getBestSNPs.R")
source("../src/hashCall.R")
source("../src/createHeatmap.R")

####################################################################
## Parameters

verbose=T

exptName="snacsExpt"
hashNames=c("CD45.26","CD45.27","CD45.28")
hashColors=c("green3","indianred2","dodgerblue3")
pvSnpThres=0.05 ## P-value for selecting best SNPs

####################################################################
## Load data
load("example.RData")

####################################################################
## Create SNACS object

snacsObj=createSNACSobject(mut=mutMat,annSNP=annSNP,annCell=annCell,exptName=exptName,hashNames=hashNames,hashColors=hashColors)
snacsObj=imputeMissingMutations(snacsObj=snacsObj,verbose=verbose)

## ------------------------------
## Cluster data to rank and select best SNPs

snacsObj=getRankedSNPs(snacsObj,col_anno_var=snacsObj$annHash$hashNames)
snacsObj=getBestSNPs(snacsObj,col_anno_var=c(snacsObj$annHash$hashNames,"clustRankedSNPs_hclust"),pvSnpThres=pvSnpThres)

## ------------------------------
## Make hash calls

## Parameters for making hash calls
prob_ecdf=0.95
propAboveBackground=0.5
probBackgndPeak=0.6
probBelowForeground=0.75

generateHashDensityPlot(snacsObj,prob_ecdf=prob_ecdf,probBackgndPeak=probBackgndPeak)
snacsObj=makeHashCall(snacsObj,prob_ecdf=prob_ecdf,probBackgndPeak=probBackgndPeak,probBelowForeground=probBelowForeground,propAboveBackground=propAboveBackground,minClustSize=500,clustComparePValue=10^-5)

## ------------------------------
## Create pretty heatmap

clustObj=createHeatmap(snacsObj,col_anno_var=c(snacsObj$annHash$hashNames,"hashCall","clustRankedSNPs_hclust"),col_anno_name=c(snacsObj$annHash$hashNames,"hash","rankedSNPsCluster"),col_dend=T,row_dend=F,outputfileName=paste0("heatmap_pvBestWithAnn_",snacsObj$exptName))

####################################################################
####################################################################
