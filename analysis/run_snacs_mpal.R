

####################################################################
####################################################################
source("../src/process_data.R")
source("../src/filter_data.R")
source("../src/create_heatmap.R")
source("../src/rank_SNPs.R")
source("../src/heatmapRelated.R")

####################################################################
## Parameters

exptName="mpal3"
hashNames=c("CD45.26","CD45.27","CD45.28")

exptName="mpal3"

dirData="../data/"
load(paste0(dirData,"m.int.",exptName,".RData"))
load(paste0(dirData,"cell_barcode_",exptName,".RData"))
datObj=process_data(mutation=mmpalDat.int,annSNP=annSNP,annCell=annCell)

#datObj=impute_data(datObj=datObj)
load(paste0(dirData,"m.int.filt.imp.",exptName,".RData"))
load(paste0(dirData,"cell_barcode_",exptName,".RData"))
datObj=list(mutation=mmpalDat.int.filt.imp,annSNP=annSNP,annCell=annCell)

datObj=rank_SNPs(datObj0,col_anno_var=hashNames,hashColors=c("green3","indianred2","dodgerblue3"))

####################################################################

####################################################################
####################################################################
