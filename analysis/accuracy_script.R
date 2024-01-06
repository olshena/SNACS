## SNACS_Manuscript
## Get accuracy estimates for each multi-sample experiment

## Generate annCell tables with information for accuracy estimates
## Generate heatmaps with truth call annotation

## Calculation of accuracy in a multi-sample experiment:
## We consider the SNPs from the multi-sample experiment that are present in all the constituent single sample experiments.

####################################################################
####################################################################

library(SNACS)
source("../analysis/functions.R")

## ------------------------------
## Impute single experiments

dirData <- "../data/"
fileList <- paste0("snacsObj_init_SNACS",c(1:4),"_unfilt.RData")
for (fId in 1:length(fileList)) {
    load(paste0(dirData,fileList[fId]))
    cat("\n\n----------------- \n")
    print(snacsObj)
    
    snacsObj <- SNACS::filterData(snacsObj=snacsObj)
    snacsObj <- SNACS::imputeMissingMutations(snacsObj=snacsObj)
    print(snacsObj)
    save(snacsObj,file=paste0(dirData,"snacsObj_imp_",snacsObj$exptName,".RData"))
}


####################################################################
####################################################################
## Run to generate initial annCell files for making truth calls

#for (snacsExpt in c("SNACS5","SNACS6","SNACS7")) {
for (snacsExpt in c("SNACS5","SNACS7")) {
#for (snacsExpt in c("SNACS6")) {
    dirData <- "../data/"
    for (numSNP in numSNPvec) {
        fName <- paste0("snacsObj_",snacsExpt,"_unfilt.RData")
        cat("\n\n----------------- Input file: ",fName,"\n",sep="")
        load(paste0(dirData,fName))
        getAnnoForAccuracy(snacsObj,hashCallMismatch=2,accuracy=F,exptNameSingleSuffix="_unfilt")
    }
}

####################################################################
####################################################################
## Make truth calls

library(SNACS)
source("../analysis/functions.R")
source("../analysis/truthCall.R")

for (exptNameSuffix in paste0("_unfilt")) {
    for (thresGenoCall in c(0.005)) {
        makeTruthCall(exptNameSuffix,thresGenoCall)
    }
}

####################################################################
####################################################################
## Run to generate heatmaps based on multiplet/singlet experiment common SNPs
## Get sens/spec/acc estimates

library(SNACS)
source("../analysis/functions.R")

dirTrueHashCall="../output/accuracy/"

#for (snacsExpt in c("SNACS5","SNACS6","SNACS7")) {
#for (snacsExpt in c("SNACS5")) {
for (snacsExpt in c("SNACS6","SNACS7")) {
dirData <- "../data/"
    fName <- paste0("snacsObj_withDoubletD_",snacsExpt,"_unfilt.RData")
    cat("\n\n----------------- Input file: ",fName,"\n",sep="")
    
    load(paste0(dirData,fName))
    fName <- paste0("_",snacsExpt,"_unfilt_",snacsExpt,"snp")
    snacsObj=getTruthCall(snacsObj,fName=fName,dirTrueHashCall=dirTrueHashCall)
    getAnnoForAccuracy(snacsObj,hashCallMismatch=2,accuracy=T,exptNameSingleSuffix="_unfilt",dirTrueHashCall=dirTrueHashCall,dirOutput=dirOutput,writeOutput="")
    createAnnotatedHeatmap(snacsObj,dirOutput=dirOutput2)
}

####################################################################
####################################################################
## Get sens/spec/acc estimates

exptNameSuffix="_unfilt"
exptNameSuffix="_hashFiltNumSNP3"

dirTrueHashCall="../output/accuracy/trueHashCall/"

for (snacsExpt in c("SNACS5","SNACS6","SNACS7")) {
    dirData <- "../data/"
    fName <- paste0("snacsObj_withDoubletD_",snacsExpt,exptNameSuffix,".RData")
    cat("\n\n----------------- ",snacsExpt,"\n",sep="")
    load(paste0(dirData,fName))
    fName <- paste0("_",snacsExpt,exptNameSuffix,"_",snacsExpt,"snp")
    snacsObj=getTruthCall(snacsObj,fName=fName,dirTrueHashCall=dirTrueHashCall)
    
    ## ----------------------------------------
    phen=getExptInfoData()
    x=table(phen$run,phen$patient)
    k=apply(x,1,function(x) {sum(x!=0)==1})
    phen=phen[which(phen$run%in%rownames(x)[k]),]

    snacsObj$annHash$exptName=phen$run[match(snacsObj$annHash$patient,phen$patient)]
    patInfo=data.frame(patTruth=c("multiplet","ambiguous",paste0("Sample.",1:nrow(snacsObj$annHash))),hash=c("Multiplet","Ambiguous",snacsObj$annHash$hashNames),pat=c("Multiplet","Ambiguous",sub("Patient ","pat",snacsObj$annHash$patient)),exptName=c("Multiplet","Ambiguous",sub("Patient ","pat",snacsObj$annHash$exptName)))
    x=snacsObj$annCell
    x$obs=x$snacsRnd2
    x$obs[which(!x$obs%in%snacsObj$annHash$hashNames)]="Multiplet"
    x$obs=patInfo$pat[match(x$obs,patInfo$hash)]
    y=x
    y$hashCall_truth[which(x$hashCall_truth%in%c("Ambiguous"))]=NA
    y1=x[which(x$hashCall_truth%in%c("Multiplet")),c("hashCall_truth","obs")]
    y0=x
    y0$hashCall_truth[which(x$hashCall_truth%in%c("Multiplet","Ambiguous"))]=NA
    print(table(x$hashCall_truth,x$obs,exclude=NULL))
    table(y1$hashCall_truth,y1$obs)
    table(y0$hashCall_truth,y0$obs)
    table(y$hashCall_truth,y$obs)
    cat("Sens/spec/acc: ",round(mean(y1$hashCall_truth==y1$obs,na.rm=T),2)," / ",round(mean(y0$hashCall_truth==y0$obs,na.rm=T),2)," / ",round(mean(y$hashCall_truth==y$obs,na.rm=T),2),"\n",sep="")
    ## ----------------------------------------
    
    #createAnnotatedHeatmap(snacsObj,dirOutput=dirOutput2)
}

####################################################################
####################################################################
