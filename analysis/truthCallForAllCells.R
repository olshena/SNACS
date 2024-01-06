## SNACS_Manuscript
## Make truth calls for all cells

####################################################################
####################################################################
library(SNACS)
library(heatmap4)
source("../analysis/functions.R")
source("../analysis/truthCall.R")

####################################################################
####################################################################
imputeMissingMutForAllCells=function(snacsObj,seed=12345,verbose=FALSE) {
    timeStamp=Sys.time()

    #if (is.na(match("filtered",snacsObj[["processLevels"]]))) stop("Run getBestSNPs() before imputing data\n")

    if (verbose) cat("\n\nImputing ",snacsObj$exptName," ...\n",sep="")
    dirData="../data/"
    
    datThis=snacsObj$mut

    missMat=matrix(F,nrow=nrow(datThis),ncol=ncol(datThis),dimnames=list(rownames(datThis),colnames(datThis))); missMat[is.na(datThis)]=T
    x=matrix(nrow=nrow(datThis),ncol=ncol(datThis),dimnames=list(rownames(datThis),colnames(datThis)))
    x[datThis==0]="F"; x[datThis==1]="T"
    cellNames=colnames(x)
    rm(datThis)
    x=t(x)
    #set.seed(seed)
    x=VIM::kNN(x,k=5,imp_var=F)
    x=as.matrix(x); x=t(x)
    datThis=matrix(nrow=nrow(x),ncol=ncol(x),dimnames=list(rownames(x),cellNames))
    datThis[x=="F"]=0
    datThis[x=="T"]=1
    rm(x)

    snacsObj[["mut"]]=datThis
    snacsObj[["missing"]]=missMat
    snacsObj[["processLevels"]]=c(snacsObj[["processLevels"]],"imputed")

    if (verbose) {
        cat("\nImputation done\n\n",sep="")
        timeStamp=c(timeStamp,Sys.time())
        print(diff(timeStamp))
    }
    
    invisible(snacsObj)
}

####################################################################
####################################################################

numSNPvec=1:5
exptIdVec=5:7

thresGenoCall=0.005

dirOutput="../output/accuracy/allCell/"
if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))

## ------------------------------
## Create single experiment annCell files

dirData="../output/accuracy/annCell/"
for (exptId in 1:4) {
    fileList=dir(dirData,pattern=paste0("SNACS",exptId))
    for (fId in 1:length(fileList)) {
        annCell=read.table(paste0(dirData,fileList[fId]),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
        fName=strsplit(fileList[fId],"_")[[1]]
        fName=paste0(paste0(fName[1:3],collapse="_"),"_allCell_",fName[4])
        write.table(annCell,file=paste0(dirOutput,"annCell/",fName),col.names=T,row.names=F, sep="\t",quote=F)
    }
}

## ------------------------------
## Make truth calls

for (numSNP in numSNPvec) {
    exptNameSuffix=paste0("_hashFiltNumSNP",numSNP)
    
    for (exptId in exptIdVec) {
        dirData <- "../data/"
        
        fName <- paste0("snacsObj_init_SNACS",exptId,exptNameSuffix,".RData")
        #fName <- paste0("snacsObj_withDoubletD_SNACS",exptId,exptNameSuffix,"_SNACS",exptId,"snp.RData")
        load(paste0(dirData,fName))
        snacsObj0=snacsObj

        fName <- paste0("snacsObj_withDoubletD_SNACS",exptId,exptNameSuffix,".RData")
        #fName <- paste0("snacsObj_withDoubletD_SNACS",exptId,exptNameSuffix,"_SNACS",exptId,"snp.RData")
        load(paste0(dirData,fName))
        snacsObj1=snacsObj
        cat("\n\n----------------- ",snacsObj$exptName,"\n",sep="")
        
        dirData="../output/accuracy/snpPosPropMat/"
        fName <- paste0("snpPosPropMat_SNACS",exptId,exptNameSuffix,"_SNACS",exptId,"snp.txt")
        snpPosPropMat=read.table(paste0(dirData,fName),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)

        dirData="../output/accuracy/annCell/"
        fName <- paste0("annCell_SNACS",exptId,exptNameSuffix,"_SNACS",exptId,"snp.txt")
        annCell=read.table(paste0(dirData,fName),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
                
        ## ------------------------------
        i=match(snacsObj1$annSNP$id,snacsObj0$annSNP$id)
        snacsObj0$mut=snacsObj0$mut[i,]
        snacsObj0$annSNP=snacsObj0$annSNP[i,]

        j=match(snacsObj0$annCell$id,snacsObj1$annCell$id); j0=which(!is.na(j)); j1=j[j0]; j00=which(is.na(j))
        
        snacsObj=snacsObj0
        snacsObj$exptName=paste0(snacsObj$exptName,"_allCell")
        snacsObj <- imputeMissingMutForAllCells(snacsObj=snacsObj)
        snacsObj$mut[,j0]=snacsObj1$mut[,j1]
        print(table(c(snacsObj$mut[,j0])==c(snacsObj1$mut[,j1]),exclude=NULL)) ## All TRUE
        print(table(c(snacsObj$mut),exclude=NULL)) ## No NA
        snacsObj0=snacsObj

        i=match(snpPosPropMat$id,snacsObj1$annSNP$id)
        snacsObj1$mut=snacsObj1$mut[i,]
        snacsObj1$annSNP=snacsObj1$annSNP[i,]

        i=match(snpPosPropMat$id,snacsObj0$annSNP$id)
        snacsObj0$mut=snacsObj0$mut[i,]
        snacsObj0$annSNP=snacsObj0$annSNP[i,]
        
        if (length(i)==1) {
            snacsObj1$mut=matrix(snacsObj1$mut,nrow=1)
            snacsObj0$mut=matrix(snacsObj0$mut,nrow=1)
        }

        ## ------------------------------
        j=match(snacsObj0$annCell$id,snacsObj1$annCell$id); j0=which(!is.na(j)); j1=j[j0]; j00=which(is.na(j))
        snacsObj=snacsObj0

        if (length(j00)!=0) {
            colId=c("id","desc")
            tmp=annCell[1:length(j00),]
            for (k in 1:ncol(tmp)) tmp[,k]=NA
            tmp[,colId]=snacsObj$annCell[j00,colId]
            genotype=as.character(snacsObj$mut[1,j00])
            if (nrow(snacsObj$mut)>1) {for (i in 2:nrow(snacsObj$mut)) {genotype=paste(genotype,snacsObj$mut[i,j00])}}
            tmp$genotype=genotype
            annCell=rbind(annCell,tmp)
        }
        
        dirOutThis=paste0(dirOutput,"annCell/")
        fName <- paste0("_",snacsObj$exptName,"_SNACS",exptId,"snp")
        if (!file.exists(dirOutThis)) dir.create(file.path(dirOutThis))
        write.table(annCell,file=paste0(dirOutput,"annCell/annCell",fName,".txt"),col.names=T,row.names=F, sep="\t",quote=F)
    }

    exptNameSuffixNew=paste0("_hashFiltNumSNP",numSNP,"_allCell")
    makeTruthCall(exptNameSuffixNew,thresGenoCall,dirData=paste0(dirOutput,"annCell/"),dirOutput=paste0(dirOutput,"trueHashCall/"))
}

## ------------------------------
## Add truth calls to annCell table

for (numSNP in numSNPvec) {
    exptNameSuffix=paste0("_hashFiltNumSNP",numSNP)
    exptNameSuffixNew=paste0("_hashFiltNumSNP",numSNP,"_allCell")
    for (exptId in exptIdVec) {
        dirData <- "../data/"
        fName <- paste0("snacsObj_withDoubletD_SNACS",exptId,exptNameSuffix,".RData")
        load(paste0(dirData,fName))
        cat("\n\n----------------- ",snacsObj$exptName,"\n",sep="")

        patInfoAll=data.frame(patTruth=c("multiplet","ambiguous",paste0("Sample.",1:nrow(snacsObj$annHash))),hash=c("Multiplet","Ambiguous",snacsObj$annHash$hashNames),pat=c("Multiplet","Ambiguous",sub("Patient ","pat",snacsObj$annHash$patient)))

        dirData="../output/accuracy/"
        fName <- paste0("annCell_SNACS",exptId,exptNameSuffix,"_SNACS",exptId,"snp.txt")
        annCell1=read.table(paste0(dirData,"annCell/",fName),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)

        fName <- paste0("annCell_SNACS",exptId,exptNameSuffixNew,"_SNACS",exptId,"snp.txt")
        annCell=read.table(paste0(dirOutput,"annCell/",fName),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
        fName <- paste0("trueHashCall_SNACS",exptId,exptNameSuffixNew,"_SNACS",exptId,"snp.txt")
        hashCall_truth=read.table(paste0(dirOutput,"trueHashCall/",fName),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
        
        annCell1=annCell
        
        annCell$hashCall_truth=hashCall_truth$hashCall_truth[match(annCell$id,hashCall_truth$id)]
        patInfo=patInfoAll[match(sort(unique(annCell$hashCall_truth)),patInfoAll$patTruth),]
        annCell$hashCall_truth=patInfo$pat[match(annCell$hashCall_truth,patInfo$patTruth)]
        
        j=match(annCell1$id,annCell$id)
        annCell2=annCell[j,]
        
        #print(table(snacs=annCell1$hashCall_truth,allCell=annCell2$hashCall_truth,exclude=NULL))
        j=which(annCell1$hashCall_truth!=annCell2$hashCall_truth | is.na(annCell1$hashCall_truth)!=is.na(annCell2$hashCall_truth))
        j=which(is.na(annCell$hashClust))
        if (length(j)!=0) {
            print(table(snacs=annCell1$hashCall_truth[j],allCell=annCell2$hashCall_truth[j],exclude=NULL))
        }
        
        fName <- paste0("annCell_SNACS",exptId,exptNameSuffixNew,"_SNACS",exptId,"snp.txt")
        if (!file.exists(dirOutThis)) dir.create(file.path(dirOutThis))
        write.table(annCell,file=paste0(dirOutput,"annCell/",fName),col.names=T,row.names=F, sep="\t",quote=F)
    }
}

####################################################################
####################################################################
