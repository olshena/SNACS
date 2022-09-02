
####################################################################
####################################################################
source("../src/heatmapRelated.R")
source("../src/process_data.R")
source("../src/filter_data.R")
source("../src/create_heatmap.R")
source("../src/rank_SNPs.R")
source("../src/getBestSNPs.R")

####################################################################
## Parameters

exptName="mpal3"
hashNames=c("CD45.26","CD45.27","CD45.28")

exptName="mpal3"

pvSnpThres=0

####################################################################
dirData="../data/"
load(paste0(dirData,"m.int.",exptName,".RData"))
load(paste0(dirData,"cell_barcode_",exptName,".RData"))
datObj=process_data(mut=mmpalDat.int,annSNP=annSNP,annCell=annCell)

#datObj=impute_data(datObj=datObj)
load(paste0(dirData,"m.int.filt.imp.",exptName,".RData"))
load(paste0(dirData,"cell_barcode_",exptName,".RData"))
datObj=list(mut=t(mmpalDat.int.filt.imp),annSNP=annSNP,annCell=annCell)

datObj=rank_SNPs(datObj,col_anno_var=hashNames,hashColors=c("green3","indianred2","dodgerblue3"))

datObj=getBestSNPs(datObj,col_anno_var=hashNames,hashColors=c("green3","indianred2","dodgerblue3"),pvSnpThres=pvSnpThres)

## ------------------------------

source("../src/hashCall.R")

prob_ecdf=0.95
propAboveBackground=0.5
probBackgndPeak=0.6
probBelowForeground=0.75

mPalIdVec=3
for (exptName in paste0("mpal",mPalIdVec)) {
    hashDensityPlot(exptName=exptName,hashNames=hashNames,prob_ecdf=prob_ecdf,probBackgndPeak=probBackgndPeak)
    hashCall(exptName=exptName,hashNames=hashNames,prob_ecdf=prob_ecdf,probBackgndPeak=probBackgndPeak,probBelowForeground=probBelowForeground,propAboveBackground=propAboveBackground,minClustSize=500,clustComparePValue=10^-5)
}
for (exptName in paste0("mpal",mPalIdVec)) {
}

####################################################################

####################################################################
####################################################################
