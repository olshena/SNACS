
argv=commandArgs(TRUE)

print(argv)
if (length(argv)!=0) {
    argv=t(sapply(argv,function(x) {
        y=strsplit(x,"=")[[1]]
        if (length(y)==1) y=c(y,"")
        y
    },USE.NAMES=F))
    
    numCell=NA
    percMult=NA
    nRep=NA
    simIdStart=NA
    
    for (k in 1:nrow(argv)) {
        switch(argv[k,1],
            "nRep"={
                nRep=as.integer(argv[k,2])
            },
            "simIdStart"={
                simIdStart=as.integer(argv[k,2])
           }
        )
    }
}


## ------------------------------
## Parameters

simIdVec=simIdStart+(1:nRep)-1
numSamVec=2:4
numCellVec=seq(2000,10000,by=2000)
propMultVec=c(0,0.05,0.1,0.15,0.20,0.25)

outputFormat="none"

fNameSuffix=""

## ------------------------------
library(SNACS)

source("run_simulation.R")
source("functions_sim.R")

dirData="../data/"
dirOutput="../output/"
dirLog="../log/"


cat("\n\n-------------------------\n")
cat("SNACS simulations: ",paste0(range(simIdVec),collapse=" - "),"\n",sep="")


## ------------------------------

runTime=Sys.time()
cat("\nStart time: ",format(runTime, "%x %X"),"\n",sep="")

## ------------------------------
# Run SNACS

## Timestamp simulation run info file with yyMMddhhmm
timeStamp=substr(gsub("-|:| ","",format(Sys.time())),3,12)

## Table to collect simulation run info
n=length(numSamVec)*length(numCellVec)*length(propMultVec)*length(simIdVec)
tmp=rep(NA,n)
runInfo=data.frame(run=rep("",n),simId=tmp,numCell=tmp,propMult=tmp,numSam=tmp,successSNACS=tmp,successDoubletD=tmp,numCellTotal=tmp,numCellSnacs=tmp,accSnacs=tmp,sensSnacs=tmp,specSnacs=tmp,correctSingle=tmp,accSnacsPlusDoubletD=tmp,sensSnacsPlusDoubletD=tmp,specSnacsPlusDoubletD=tmp,correctSingleSnacsPlusDD=tmp,hashBgnd=rep("",n))

kRun=1
#kRun=73

for (simId in simIdVec) {
    for (numCell in numCellVec) {
        for (propMult in propMultVec) {
            for (numSam in numSamVec) {
                runInfo$run[kRun]=kRun
                runInfo$simId[kRun]=simId
                runInfo$numCell[kRun]=numCell
                runInfo$propMult[kRun]=propMult
                runInfo$numSam[kRun]=numSam

                exptName=paste0("sim",simId,"_",numSam,"samples_",numCell,"cells_",propMult,"propMult",fNameSuffix)
                load(file=paste0(dirData,"snacsObj_init_",exptName,".RData"))
                runInfo$numCellTotal[kRun]=nrow(snacsObj$annCell)
                ## ------------------------------
                cat("\n\n----------------- ",snacsObj$exptName,"\n",sep="")

                ## Get true call for all cells in simulated data
                callObsAll=sapply(snacsObj$annCell$truth,modifyCallName,USE.NAMES=F)
                names(callObsAll)=snacsObj$annCell$id
                
                ## Run SNACS
                snacsObj <- try(SNACS::runSNACS(snacsObj=snacsObj,hashThreshold=NA,outputFormat=outputFormat))

                if ("SNACSList"%in%class(snacsObj)) {

                    snacsObj$missing=snacsObj$centroid=snacsObj$dist2centroidMat=NULL
                    
                    save(snacsObj,file=paste0(dirData,"snacsObj_",snacsObj$exptName,".RData"))

                    runInfo$numCellSnacs[kRun]=nrow(snacsObj$annCell)
                    runInfo$hashBgnd[kRun]=paste(snacsObj$annHashBackgnd$method,collapse=",")

                    ## Consider true call for cells in final snacsObj only
                    callObs=callObsAll[match(snacsObj$annCell$id,names(callObsAll))]

                    ## Get SNACS call for the cells in simulated data
                    callPred=sapply(snacsObj$annCell$snacsRnd2,modifyCallName,USE.NAMES=F)
                    
                    ## Get acc/sens/spec estimates for SNACS calls
                    runInfo$accSnacs[kRun]=mean(callPred==callObs,na.rm=T)
                    j=which(callObs=="multiplet")
                    runInfo$sensSnacs[kRun]=sum(callPred[j]==callObs[j],na.rm=T)/length(j)
                    j=which(callObs%in%snacsObj$annHash$hashNames)
                    runInfo$specSnacs[kRun]=sum(callPred[j]==callObs[j],na.rm=T)/length(j)
                    j=which(callObs%in%snacsObj$annHash$hashNames & callPred%in%snacsObj$annHash$hashNames)
                    runInfo$correctSingle[kRun]=sum(callPred[j]==callObs[j],na.rm=T)/length(j)
                    
                    runInfo$successSNACS[kRun]=1

                    ## Run SNACS+doubletD
                    snacsObj <- try(SNACS::runSNACSplusDoubletD(snacsObj))
                    if ("SNACSList"%in%class(snacsObj)) {
                        save(snacsObj,file=paste0(dirData,"snacsObj_withDoubletD_",snacsObj$exptName,".RData"))
                        
                        res=getColVarInfo(snacsObj$annCell,snacsObj$annHash)
                        if (outputFormat!="none") {
                            clustObj=createHeatmapForSim(snacsObj,colVarInfo=res$col_var_info,cell_anno_var=res$cellAnnoVar,cell_anno_name=res$cellAnnoName,col_dend=TRUE,row_dend=FALSE,outputFileName=paste0("heatmap_",snacsObj$exptName),outputFormat="pdf")
                            #clustObj=SNACS::createHeatmap(snacsObj,cell_anno_var=c("snacsRnd2","snacsRnd1","clustBestSNPs_hclust",snacsObj$annHash$hashNames),cell_anno_name=c("snacsRnd2","snacsRnd1","bestSNPsCluster",snacsObj$annHash$hashNames),col_dend=TRUE,row_dend=FALSE,outputFileName=paste0("heatmap_",snacsObj$exptName),outputFormat="pdf")
                        }

                        ## ----------------------------------------------
                        ## Get SNACS+doubletD call for the cells in simulated data
                        callPred=sapply(snacsObj$annCell$snacsPlusDoubletD,modifyCallName,USE.NAMES=F)
                        
                        ## Get acc/sens/spec estimates for SNACS+doubletD calls
                        runInfo$accSnacsPlusDoubletD[kRun]=mean(callPred==callObs,na.rm=T)
                        j=which(callObs=="multiplet")
                        runInfo$sensSnacsPlusDoubletD[kRun]=sum(callPred[j]==callObs[j],na.rm=T)/length(j)
                        j=which(callObs%in%snacsObj$annHash$hashNames)
                        runInfo$specSnacsPlusDoubletD[kRun]=sum(callPred[j]==callObs[j],na.rm=T)/length(j)
                        j=which(callObs%in%snacsObj$annHash$hashNames & callPred%in%snacsObj$annHash$hashNames)
                        runInfo$correctSingleSnacsPlusDD[kRun]=sum(callPred[j]==callObs[j],na.rm=T)/length(j)

                        runInfo$successDoubletD[kRun]=1
                    } else {
                        res=getColVarInfo(snacsObj$annCell,snacsObj$annHash)
                        if (outputFormat!="none") {
                            clustObj=createHeatmapForSim(snacsObj,colVarInfo=res$col_var_info,cell_anno_var=res$cellAnnoVar,cell_anno_name=res$cellAnnoName,col_dend=TRUE,row_dend=FALSE,outputFileName=paste0("heatmap_",snacsObj$exptName),outputFormat="pdf")
                            #clustObj=SNACS::createHeatmap(snacsObj,cell_anno_var=c("snacsRnd2","snacsRnd1","clustBestSNPs_hclust",snacsObj$annHash$hashNames),cell_anno_name=c("snacsRnd2","snacsRnd1","bestSNPsCluster",snacsObj$annHash$hashNames),col_dend=TRUE,row_dend=FALSE,outputFileName=paste0("heatmap_",snacsObj$exptName),outputFormat="pdf")
                        }
                        
                        runInfo$successDoubletD[kRun]=0
                    }
                } else {
                    runInfo$successSNACS[kRun]=0
                }
                kRun=kRun+1
            }
        }
    }
}

if (T) {
    ## Write runInfo table to file
    tbl=runInfo
    for (k in c("accSnacs","sensSnacs","specSnacs","accSnacsPlusDoubletD","sensSnacsPlusDoubletD","specSnacsPlusDoubletD","correctSingle","correctSingleSnacsPlusDD")) tbl[,k]=round(tbl[,k],3)
    write.table(tbl,paste0(dirOutput,"runInfo_SNACS_",timeStamp,".txt"),sep="\t",col.names=T,row.names=F,quote=F)
}
#tbl

runTime=c(runTime,Sys.time())
cat("End time: ",format(runTime[2], "%x %X"),"\n",sep="")
print(diff(runTime))

## ------------------------------

####################################################################
####################################################################
