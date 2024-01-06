## SNACS_filter
## Get accuracy estimates for each multi-sample experiment

## Calculation of accuracy in a multi-sample experiment:
## Generate heatmaps with truth call annotation

####################################################################
####################################################################


library(SNACS)
library(heatmap4)
source("../analysis/functions.R")
source("../analysis/truthCall.R")

####################################################################
####################################################################

## ------------------------------
## Impute single experiments

numSNPvec=1:5

for (exptName in paste0("SNACS",1:4)) {
    ## -------------------------------------
    dirData <- "../data/"
    outputFormat="pdf"

    load(file=paste0(dirData,"snacsObj_init_",exptName,"_unfilt.RData"))
    snacsObj
    for (numSNP in numSNPvec) {
        snacsObj$exptName=paste0(exptName,"_hashFiltNumSNP",numSNP)
        save(snacsObj,file=paste0(dirData,"snacsObj_init_",snacsObj$exptName,".RData"))
    }

    ####################################################################
    ####################################################################

    snacsObj$exptName=paste0(exptName,"_hashFilt")

    snacsObj <- SNACS::filterData(snacsObj=snacsObj,proportionMutatedPerSNP=c(0.05,0.95))

    snacsObj <- SNACS::imputeMissingMutations(snacsObj=snacsObj)

    for (numSNP in numSNPvec) {
        snacsObj$exptName=paste0(exptName,"_hashFiltNumSNP",numSNP)
        save(snacsObj,file=paste0(dirData,"snacsObj_imp_",snacsObj$exptName,".RData"))
    }
}

####################################################################
####################################################################
## Run to generate initial annCell files for making truth calls

numSNPvec=1:5

for (numSNP in numSNPvec) {
    for (snacsExpt in c("SNACS5","SNACS6","SNACS7")) {
        dirData <- "../data/"
        fName <- paste0("snacsObj_",snacsExpt,"_hashFiltNumSNP",numSNP,".RData")
        cat("\n\n----------------- Input file: ",fName,"\n",sep="")
        load(paste0(dirData,fName))
        getAnnoForAccuracy(snacsObj,hashCallMismatch=2,accuracy=T,exptNameSingleSuffix=paste0("_hashFiltNumSNP",numSNP))
    }
}

####################################################################
####################################################################
## Make truth calls

numSNPvec=1:5
thresGenoCall=0.005

for (numSNP in numSNPvec) {
    exptNameSuffix=paste0("_hashFiltNumSNP",numSNP)
    makeTruthCall(exptNameSuffix,thresGenoCall)
}

####################################################################
####################################################################
## Run to generate heatmaps based on multiplet/singlet experiment common SNPs

numSNPvec=1:5
thresGenoCall=0.005

dirTrueHashCall=paste0("../output/accuracy/trueHashCall/")
dirOutput=paste0("../output/accuracy/")
dirOutput2=paste0("../output/accuracy/heatmapAllSNP/")
if (!file.exists(dirOutput2)) dir.create(file.path(dirOutput2))
for (numSNP in numSNPvec) {
    for (snacsExpt in c("SNACS5","SNACS6","SNACS7")) {
    #for (snacsExpt in c("SNACS5")) {
        dirData <- "../data/"
        fName <- paste0("snacsObj_withDoubletD_",snacsExpt,"_hashFiltNumSNP",numSNP,".RData")
        cat("\n\n----------------- Input file: ",fName,"\n",sep="")
        
        if (snacsExpt=="SNACS7" & numSNP==5 & thresGenoCall==0.002) next
        
        load(paste0(dirData,fName))
        fName <- paste0("_",snacsExpt,"_hashFiltNumSNP",numSNP,"_",snacsExpt,"snp")
        snacsObj=getTruthCall(snacsObj,fName=fName,dirTrueHashCall=dirTrueHashCall)
        getAnnoForAccuracy(snacsObj,hashCallMismatch=2,accuracy=T,exptNameSingleSuffix=paste0("_hashFiltNumSNP",numSNP),dirTrueHashCall=dirTrueHashCall,dirOutput=dirOutput,writeOutput="")
        createAnnotatedHeatmap(snacsObj,dirOutput=dirOutput2)
    }
}

####################################################################
####################################################################
## Get sens/spec/acc estimates

numSNPvec=1:5
thresGenoCall=0.005

dirTrueHashCall=paste0("../output/accuracy/trueHashCall/")
dirOutput=paste0("../output/accuracy/")
dirOutput2=paste0("../output/accuracy/heatmapAllSNP/")
if (!file.exists(dirOutput2)) dir.create(file.path(dirOutput2))
for (numSNP in numSNPvec) {
    #for (snacsExpt in c("SNACS5")) {
    for (snacsExpt in c("SNACS5","SNACS6","SNACS7")) {
        dirData <- "../data/"
        fName <- paste0("snacsObj_withDoubletD_",snacsExpt,"_hashFiltNumSNP",numSNP,".RData")
        cat("\n\n----------------- ",snacsExpt,", ",numSNP," SNP-pairs\n",sep="")
        load(paste0(dirData,fName))
        fName <- paste0("_",snacsExpt,"_hashFiltNumSNP",numSNP,"_",snacsExpt,"snp")
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
        
        createAnnotatedHeatmap(snacsObj,dirOutput=dirOutput2)
    }
}

####################################################################
####################################################################
