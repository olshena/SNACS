
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
## Get simulated SNACS objects

load("../data/data2SimulateFrom.RData")

for (simId in simIdVec) {
    for (numCell in numCellVec) {
        for (propMult in propMultVec) {
            for (numSam in numSamVec) {
                set.seed(simId)
                switch(as.character(numSam),
                    "2"={exptSingleId=1:2
                        fooThis=sim(datasets.list=list(s1,s2),n.cells=numCell,prop.doublets=propMult,improper.multiplier=0.05)},
                    "3"={
                        exptSingleId=2:4
                        fooThis=sim(datasets.list=list(s2,s3,s4),n.cells=numCell,prop.doublets=propMult,improper.multiplier=0.05)
                    },
                    "4"={exptSingleId=1:4
                        fooThis=sim(datasets.list=list(s1,s2,s3,s4),n.cells=numCell,prop.doublets=propMult,improper.multiplier=0.05)}
                )

                exptName=paste0("sim",simId,"_",numSam,"samples_",numCell,"cells_",propMult,"propMult",fNameSuffix)
                hashColors <- c("indianred2","green3","dodgerblue3","magenta3")
                rownames(fooThis$hash.output)=paste0("TS.",exptSingleId)
                if (numSam==3) {
                    truth=truthName2snacsCallName(fooThis$truth,hashNameInit=as.character(1:3),hashNameFinal=rownames(fooThis$hash.output))
                } else {
                    truth=truthName2snacsCallName(fooThis$truth)
                }
                snacsObj=SNACSList(mut=fooThis$final.genotype.output,hashes=fooThis$hash.output,exptName=exptName,hashColors=hashColors[exptSingleId],
                                      depthTotal=fooThis$final.depthtotal.output,depthAlt=fooThis$final.depthalt.output,annCell=data.frame(desc=paste0("c",1:ncol(fooThis$final.genotype.output)),truth=truth),annSNP=data.frame(desc=paste0("s",1:nrow(fooThis$final.genotype.output))))
                rm(fooThis,truth)
                
                save(snacsObj,file=paste0(dirData,"snacsObj_init_",snacsObj$exptName,".RData"))
                
                rm(snacsObj)

            }
            
        }
    }
}

## ------------------------------
runTime=c(runTime,Sys.time())
cat("\n\nEnd time: ",format(runTime[2], "%x %X"),"\n",sep="")
print(diff(runTime))

## ------------------------------

####################################################################
####################################################################
