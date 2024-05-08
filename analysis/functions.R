
####################################################################
####################################################################

## ------------------------------
## Get spherical k-means based distance object
getSKmeansDist <- function(x) {stats::as.dist(skmeans::skmeans_xdist(x))}

## ------------------------------
## Set experiment name
setExptName <- function(snacsObj,exptName) {
    snacsObj$exptName=exptName
    invisible(snacsObj)
}

## ------------------------------
## Get experiment infomation from SNACS metadata Excel file
getExptInfoData <- function(fileName="SNACS_Metadata.csv") {
phen = read.csv(fileName)
phen=as.data.frame(phen,stringsAsFactors=F)
names(phen)[match(c("Experiment","Patient","Hash"),names(phen))]=
        c("run","patient","hash")
phen$hash=sub("-",".",phen$hash)
phen$patId=as.integer(as.factor(phen$patient))
phen$hashId=as.integer(as.factor(phen$hash))

    invisible(phen)
}

## ------------------------------
## Add hash calls from HTOdemux to SNACS object
getHTOdemuxCall <- function(snacsObj,dirName="../data/SNACS_HTOdemux/") {
    #fName=paste0(sub("_unfilt","",snacsObj$exptName),"_HTOdemux.csv")
    fName=paste0(strsplit(snacsObj$exptName,"_")[[1]][1],"_HTOdemux.csv")
    dat=read.table(paste0(dirName,fName),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
    names(dat)=c("id","hashCall")
    for (k in 1:ncol(dat)) dat[,k]=gsub("\"","",dat[,k])
    dat$hashCall=gsub("-",".",dat$hashCall)
    snacsObj$annCell$desc=gsub("-",".",snacsObj$annCell$desc)
    j=match(snacsObj$annCell$desc,dat$id); j1=which(!is.na(j)); j2=j[j1]
    snacsObj$annCell$HTOdemux="Missing"
    snacsObj$annCell$HTOdemux[j1]=dat$hashCall[j2]
    #snacsObj$annCell$HTOdemux[which(snacsObj$annCell$HTOdemux=="")]=NA
    
    invisible(snacsObj)
}

## ------------------------------
## Create SNACSList objects
createInitialSNACSobject <- function(fileName,dirName="../data/",hashColors,dirDepth=NULL,verbose=T) {
    cat("\n\n----------------- ",fileName,"\n")
    dat=read.table(paste0(dirName,fileName),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
    
    j=grep("TS.",colnames(dat),fixed=T)
    j=c(j[1]-1,j)
    if (verbose) {
        cat("Dimension of data file\n")
        print(dim(dat))
        cat("Last column containing SNP data and columns containing hash data\n")
        print(cbind(column=j,name=colnames(dat)[j]))
    }
    
    k=grep("TS.",colnames(dat),fixed=T)
    hashMat=t(as.matrix(dat[,k]))
    annCell=data.frame(id=paste0("cell",1:nrow(dat)),desc=dat[,1],stringsAsFactors=F)
    i=c(1,grep("TS.",colnames(dat),fixed=T))
    annSNP=data.frame(id=paste0("snp",1:(ncol(dat)-length(i))),desc=names(dat)[-i],stringsAsFactors=F)
    mutMat=dat[,-i]
    rm(dat)
    mutMat=as.matrix(mutMat)
    mutMat=t(mutMat)
    mutMat[mutMat==3]=NA
    colnames(mutMat)=colnames(hashMat)=annCell$id
    rownames(mutMat)=annSNP$id
    nm=read.table(paste0(dirName,fileName),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=1)
    nm=unlist(nm[1,])
    annSNP$desc=nm[-i]
    rm(nm)

    phen=getExptInfoData()
    
    exptName <- sub("_redo","",sub("ered.csv","",fileName))
    samInfo <- phen[which(phen$run==strsplit(exptName,"_")[[1]][1]),]
    samInfo <- samInfo[match(rownames(hashMat),samInfo$hash),]
    hashColorsThis <- hashColors[samInfo$hashId]
    annHash <- data.frame(patient=samInfo$patient,patColors=hashColors[samInfo$patId],stringsAsFactors=F)

    if (!is.null(dirDepth)) {
        snpId1=sapply(annSNP$desc,function(x) {
            y=gsub("*",".",gsub(":|/|<|>|-",".",x),fixed=T)
            if (substr(y,1,1)%in%".") y=paste0("X",y)
            y
        },USE.NAMES=F)
        cellId1=gsub(":|/",".",annCell$desc)

        df_total <- read.csv(paste0(dirDepth,substr(exptName,1,6),"_total_depth.csv"))
        df_alt <- read.csv(paste0(dirDepth,substr(exptName,1,6),"_alt_depth.csv"))
        if (any(rownames(df_total)!=rownames(df_alt))) {stop("\nrownames(df_total)!=rownames(df_alt)\n")}
        if (any(colnames(df_total)!=colnames(df_alt))) {stop("\ncolnames(df_total)!=colnames(df_alt)\n")}
        snpId2=colnames(df_total)[-1]
        cellId2=df_total[,1]
        
        if (any(snpId1!=snpId2)) {stop("\nsnpId1!=snpId2\n")}
        if (any(cellId1!=cellId2)) {stop("\ncellId1!=cellId2\n")}

        depthTotal=df_total[,-1]
        depthAlt=df_alt[,-1]
        depthTotal=as.matrix(depthTotal)
        depthAlt=as.matrix(depthAlt)
        depthTotal=t(depthTotal)
        depthAlt=t(depthAlt)
        rownames(depthTotal)=rownames(depthAlt)=rownames(mutMat)
        colnames(depthTotal)=colnames(depthAlt)=colnames(mutMat)
    } else {
        depthTotal=depthAlt=NULL
    }
    snacsObj <- SNACS::SNACSList(mut=mutMat,hashes=hashMat,exptName=exptName,depthTotal=depthTotal,depthAlt=depthAlt,hashColors=hashColorsThis,annCell=annCell,annSNP=annSNP,annHash=annHash)
    snacsObj
    dirData <- "../data/"
    fName <- sub("_redo","",exptName)
    save(snacsObj,file=paste0(dirData,"snacsObj_init_",fName,".RData"))
    
    invisible(snacsObj)
}

## ------------------------------
## Set observations with low quality and depth to missing
setLowQuality2Missing <- function(snacsObj,thres=10) {
    exptName=strsplit(snacsObj$exptName,"_")[[1]][1]
    
    datQ=read.table(paste0("../data/Quality and Depth/",exptName,"_quality.csv"),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
    datD=read.table(paste0("../data/Quality and Depth/",exptName,"_total_depth.csv"),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
    descSNP=read.table(paste0("../data/Quality and Depth/",exptName,"_quality.csv"),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=1)
    descSNP=unlist(descSNP)[-1]

    i=match(snacsObj$annSNP$desc,descSNP)+1
    j=match(snacsObj$annCell$desc,datD[,1])
    datQ=datQ[j,i]
    datD=datD[j,i]
    
    snacsObj$mut[datQ<thres]=NA
    snacsObj$mut[datD<thres]=NA
    
    invisible(snacsObj)
}

## ------------------------------
## Create low quality and depth matrices for an experiment
getQualityMatrixForExpt <- function(snacsObj) {
    exptName=strsplit(snacsObj$exptName,"_")[[1]][1]
    
    datQ=read.table(paste0("../data/Quality and Depth/",exptName,"_quality.csv"),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
    datD=read.table(paste0("../data/Quality and Depth/",exptName,"_total_depth.csv"),sep=",",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
    descSNP=read.table(paste0("../data/Quality and Depth/",exptName,"_quality.csv"),sep=",",h=F,quote="",comment.char="",as.is=T,fill=T,nrow=1)
    descSNP=unlist(descSNP)[-1]

    i=match(snacsObj$annSNP$desc,descSNP)+1
    j=match(snacsObj$annCell$desc,datD[,1])
    datQ=datQ[j,i]
    datD=datD[j,i]
    rownames(datQ)=rownames(datD)=snacsObj$annCell$id
    colnames(datQ)=colnames(datD)=snacsObj$annSNP$id
    
    datQ=as.matrix(datQ)
    datD=as.matrix(datD)

    datQ=t(datQ)
    datD=t(datD)

    datObj=list(datQ=datQ,datD=datD)
    invisible(datObj)
}

## ------------------------------
## Generate heatmap with annotation
plotHashPair=function(snacsObj,hashCallMismatch=3,dirData2="../output/accuracy/",dirOutput="../output/accuracy/") {
    snacsExpt=strsplit(snacsObj$exptName,"_")[[1]][1]
    fName=paste0("_",snacsObj$exptName,"_",snacsExpt,"snp")
    annCell=read.table(paste0(dirData2,"annCell",fName,".txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
    annCell=annCell[match(snacsObj$annCell$id,annCell$id),]

    colIdMismatch=paste0("snacsRnd",hashCallMismatch,"_truth_mismatch")
    
    ## Plot hash pairs
    pdf(paste0(dirOutput,"plotHashPair_",snacsObj$exptName,"_hashRnd",hashCallMismatch,".pdf"))
    #par(mfrow=c(3,3))
    for (hId1 in 1:(nrow(snacsObj$annHash)-1)) {
        hashId1=snacsObj$annHash$hashNames[hId1]
        for (hId2 in (hId1+1):nrow(snacsObj$annHash)) {
            hashId2=snacsObj$annHash$hashNames[hId2]
            
            header=snacsObj$exptName
            if (hashCallMismatch<3) header=paste0(header,", hashRnd3 - x")
            header=paste0(header,", hashRnd",hashCallMismatch," - \nred - mult called ",hashId1,", cyan - mult ",hashId2,", yellow - mult mult")

            plot(annCell[,hashId1],annCell[,hashId2],pch=16,col="grey50",cex=0.75,main=header,xlab=hashId1,ylab=hashId2)

            patId=sub("Patient ","pat",snacsObj$annHash$patient[which(snacsObj$annHash$hashNames==hashId1)])
            j=which(annCell[,colIdMismatch]==patId)
            points(annCell[j,hashId1],annCell[j,hashId2],pch=16,col="red",cex=0.75)
            patId=sub("Patient ","pat",snacsObj$annHash$patient[which(snacsObj$annHash$hashNames==hashId2)])
            j=which(annCell[,colIdMismatch]==patId)
            points(annCell[j,hashId1],annCell[j,hashId2],pch=16,col="cyan",cex=0.75)
            j=which(annCell[,colIdMismatch]=="Multiplet")
            points(annCell[j,hashId1],annCell[j,hashId2],pch=16,col="yellow",cex=0.75)

            ## Mark the points called a multiplet in round 3 but not in round 2
            if (hashCallMismatch==2) {
                j=which(annCell$snacsRnd2_truth_mismatch!="Multiplet" & annCell$snacsRnd3_truth_mismatch=="Multiplet")
                points(annCell[j,hashId1],annCell[j,hashId2],pch="x",col="black",cex=0.75)
            }
        }
    }
    dev.off()
    #print(table(rnd3=annCell$snacsRnd3_truth_mismatch,rnd2=annCell$snacsRnd2_truth_mismatch))
}

## ------------------------------
## Generate heatmap with annotation
createAnnotatedHeatmap <- function(snacsObj,dirOutput="../output/heatmap/") {
    snacsObj$annCell$patient=ifelse(sum(!duplicated(snacsObj$annHash$patient))==1,snacsObj$annHash$patient[1],"Doublet")
    j=match(snacsObj$annCell$snacsRnd1,snacsObj$annHash$hashNames); j1=which(!is.na(j)); j2=j[j1]
    snacsObj$annCell$patient[j1]=snacsObj$annHash$patient[j2]
    snacsObj$annCell$patient=gsub("Patient ","pat",snacsObj$annCell$patient)
    
    snacsRndId=1:2
    
    col_var_info=list()
    hashInfo=snacsObj$annHash
    hashInfo$patient=gsub("Patient ","pat",hashInfo$patient)
    grpUniq1=sort(unique(snacsObj$annCell$snacsRnd1))
    grpUniq1=grpUniq1[which(!grpUniq1%in%hashInfo$hashNames)]
    grpUniq2=sort(unique(snacsObj$annCell$snacsRnd2))
    grpUniq2=grpUniq2[which(!grpUniq2%in%hashInfo$hashNames)]
    grpUniqSD=c()
    if ("doubletD"%in%names(snacsObj$annCell)) {
        #col_var_info[["doubletD"]]=list(color=c("navy","gray85"),level=c("Singlet","Doublet"))
        col_var_info[["doubletD"]]=list(color=c("navy","white"),level=c("Singlet","Doublet"))
    }
    if ("snacsPlusDoubletD"%in%names(snacsObj$annCell)) {
        grpUniqSD=sort(unique(snacsObj$annCell$snacsPlusDoubletD))
        grpUniqSD=grpUniqSD[which(!grpUniqSD%in%hashInfo$hashNames)]
    }
    grpUniqHTOd=c()
    if ("HTOdemux"%in%names(snacsObj$annCell)) {
        grpUniqHTOd=sort(unique(snacsObj$annCell$HTOdemux))
        grpUniqHTOd=grpUniqHTOd[which(!grpUniqHTOd%in%c(hashInfo$hashNames,"Negative","Missing"))]
    }
    grpUniq=sort(unique(c(grpUniq1,grpUniq2,grpUniqSD,grpUniqHTOd)))
    if (length(grpUniq)!=0) {
        colVec=gray(1:length(grpUniq)/length(grpUniq))
        colVec=colVec[!colVec%in%hashInfo$hashColors]
        colVec=colVec[1:length(grpUniq)]
        k=which(grpUniq%in%grpUniq1)
        col_var_info[["snacsRnd1"]]=list(color=c(snacsObj$annHash$hashColors,colVec[k]),level=c(snacsObj$annHash$hashNames,grpUniq[k]))
        k=which(grpUniq%in%grpUniq2)
        col_var_info[["snacsRnd2"]]=list(color=c(snacsObj$annHash$hashColors,colVec[k]),level=c(snacsObj$annHash$hashNames,grpUniq[k]))
        if ("snacsPlusDoubletD"%in%names(snacsObj$annCell)) {
            k=which(grpUniq%in%grpUniqSD)
            col_var_info[["snacsPlusDoubletD"]]=list(color=c(snacsObj$annHash$hashColors,colVec[k]),level=c(snacsObj$annHash$hashNames,grpUniq[k]))
        }
        if ("HTOdemux"%in%names(snacsObj$annCell)) {
            k=which(grpUniq%in%grpUniqHTOd)
            col_var_info[["hashHTOdemux"]]=list(color=c(snacsObj$annHash$hashColors,colVec[k],"gold","brown"),level=c(snacsObj$annHash$hashNames,grpUniq[k],"Negative","Missing"))
            k=which(col_var_info[["hashHTOdemux"]]$level%in%snacsObj$annCell$HTOdemux)
            col_var_info[["hashHTOdemux"]]=list(color=col_var_info[["hashHTOdemux"]]$color[k],level=col_var_info[["hashHTOdemux"]]$level[k])
        }
        grpUniqP=unique(snacsObj$annCell$patient)
        grpUniqP=grpUniqP[which(!grpUniqP%in%hashInfo$patient)]
        if (length(grpUniqP)!=0) {
            colVec=rainbow(length(grpUniqP)+nrow(hashInfo))
            colVec=gray(1:length(grpUniqP)/length(grpUniqP))
            colVec=colVec[!colVec%in%hashInfo$patColors]
            colVec=colVec[1:length(grpUniqP)]
            col_var_info[["patient"]]=list(color=c(snacsObj$annHash$patColors,colVec),level=c(hashInfo$patient,grpUniqP))
        }
    } else {
        col_var_info[["snacsRnd1"]]=list(color=c(snacsObj$annHash$hashColors,"cyan2"),level=c(snacsObj$annHash$hashNames,"Doublet"))
        col_var_info[["patient"]]=list(color=c(snacsObj$annHash$patColors,"cyan2"),level=c(gsub("Patient ","pat",snacsObj$annHash$patient),"Doublet"))
    }
    col_var_info=getColVarInfo(snacsObj,snacsRndId)
    
    nameRow=rep("",nrow(snacsObj$annSNP))
    nameRow=snacsObj$annSNP$desc

    col_dend=T; row_dend=F
    zlm=c(0,1)
    #heatmap_color=c("orangered1","royalblue3","white")
    heatmap_color=c("navy","gray85","orangered1")
    input_legend=T
    #annCell=cbind(hash=snacsObj$annCell$snacsRnd1,clustSplitRnd1=snacsObj$annCell$clustSplitRnd1,subClustRnd1=snacsObj$annCell$subClustRnd1,as.data.frame(t(snacsObj$hashes)),stringsAsFactors=F)
    if (length(grpUniq)!=0) {
        #tbl=100*snacsObj$dist2centroidMat
        tbl=snacsObj$dist2centroidMat
        tbl=as.data.frame(t(tbl))
        colnames(tbl)=paste0("dist2centr_",colnames(tbl))
        annCell=cbind(snacsRnd1=snacsObj$annCell$snacsRnd1,clustSplitRnd1=snacsObj$annCell$clustSplitRnd1,snacsRnd2=snacsObj$annCell$snacsRnd2,clustSplitRnd2=snacsObj$annCell$clustSplitRnd2,as.data.frame(t(snacsObj$hashes)),stringsAsFactors=F)
        #annCell=cbind(hash=snacsObj$annCell$snacsRnd1,clustSplitRnd1=snacsObj$annCell$clustSplitRnd1,snacsRnd2=snacsObj$annCell$snacsRnd2,clustSplitRnd2=snacsObj$annCell$clustSplitRnd2,as.data.frame(t(snacsObj$hashes)),tbl,stringsAsFactors=F)
        if ("HTOdemux"%in%names(snacsObj$annCell)) {
            annCell=cbind(hashHTOdemux=snacsObj$annCell$HTOdemux,annCell)
        }
        if ("doubletD"%in%names(snacsObj$annCell)) {
            annCell=cbind(doubletD=snacsObj$annCell$doubletD,annCell)
        }
        if ("snacsPlusDoubletD"%in%names(snacsObj$annCell)) {
            annCell=cbind(snacsPlusDoubletD=snacsObj$annCell$snacsPlusDoubletD,annCell)
        }
        if ("hashCall_truth"%in%names(snacsObj$annCell)) {
            annCell=cbind(truth=snacsObj$annCell$hashCall_truth,annCell)
        }
        nm=c("truth","hashHTOdemux","doubletD","snacsPlusDoubletD","snacsRnd2","clustSplitRnd2","snacsRnd1","clustSplitRnd1",rev(snacsObj$annHash$hashNames))
        nm=nm[which(nm%in%names(annCell))]
        annCell=annCell[,nm]

        if (F) {
            clustObjThis=snacsObj$hclustObj_bestSNPs
            j=clustObjThis$order
            tbl=snacsObj$annCell[j,]
            x=tbl$subClustRnd2
            y=tbl$clustSplitRnd2
            j2=which(x==x[length(x)])
            j1=which(x==x[min(j2)-1])
            j1=c(j1,which(x==x[min(j1)-1]))
            annCell$snacsRnd2=""
            annCell$snacsRnd2[j[j1]]="Doublet"
        }
    } else {
        annCell=cbind(snacsRnd1=snacsObj$annCell$snacsRnd1,clustSplitRnd1=snacsObj$annCell$clustSplitRnd1,as.data.frame(t(snacsObj$hashes[nrow(snacsObj$hashes):1,])),stringsAsFactors=F)
    }
    for (k in which(names(annCell)%in%c("truth"))) annCell[is.na(annCell[,k]),k]=""
    #annCell=cbind(hash=snacsObj$annCell$snacsRnd1,clustSplitRnd1=snacsObj$annCell$clustSplitRnd1,clustSplitRnd2=snacsObj$annCell$clustSplitRnd2,as.data.frame(t(snacsObj$hashes)),stringsAsFactors=F)
    #png(paste0("heatmap_",snacsObj$exptName,".png"))
    pdf(paste0(dirOutput,"heatmap_",snacsObj$exptName,".pdf"))
    x2=heatmap4::generate_heatmap(x=snacsObj$mut,distfun=distfun,col_clust=snacsObj$hclustObj_bestSNPs,col_dend=T,row_dend=row_dend,col_info=annCell,col_anno=T,zlm=zlm,heatmap_color=heatmap_color,col_var_info=col_var_info,h_title=snacsObj$exptName,row_lab=T,row_lab_vtr=nameRow,input_legend=input_legend)
    dev.off()
}

####################################################################
####################################################################

getColVarInfo=function(snacsObj,snacsRndId) {
    patInfo=cbind(c(gsub("Patient ","pat",snacsObj$annHash$patient),"Multiplet","Ambiguous",""),c(snacsObj$annHash$patColors,"grey50","white","white"))

    col_var_info=list()
    hashInfo=snacsObj$annHash
    hashInfo$patient=gsub("Patient ","pat",hashInfo$patient)
    grpUniq1=sort(unique(snacsObj$annCell$snacsRnd1))
    grpUniq1=grpUniq1[which(!grpUniq1%in%hashInfo$hashNames)]
    grpUniq2=sort(unique(snacsObj$annCell$snacsRnd2))
    grpUniq2=grpUniq2[which(!grpUniq2%in%hashInfo$hashNames)]
    grpUniq3=c()
    if ("snacsRnd3"%in%names(snacsObj$annCell)) {
        grpUniq3=sort(unique(snacsObj$annCell$snacsRnd3))
        grpUniq3=grpUniq3[which(!grpUniq3%in%hashInfo$hashNames)]
    }
    grpUniqSD=c()
    if ("doubletD"%in%names(snacsObj$annCell)) {
        #col_var_info[["doubletD"]]=list(color=c("navy","gray85"),level=c("Singlet","Doublet"))
        col_var_info[["doubletD"]]=list(color=c("navy","white"),level=c("Singlet","Doublet"))
    }
    if ("snacsPlusDoubletD"%in%names(snacsObj$annCell)) {
        grpUniqSD=sort(unique(snacsObj$annCell$snacsPlusDoubletD))
        grpUniqSD=grpUniqSD[which(!grpUniqSD%in%hashInfo$hashNames)]
    }
    grpUniqHTOd=c()
    if ("HTOdemux"%in%names(snacsObj$annCell)) {
        grpUniqHTOd=sort(unique(snacsObj$annCell$HTOdemux))
        grpUniqHTOd=grpUniqHTOd[which(!grpUniqHTOd%in%c(hashInfo$hashNames,"Negative","Missing"))]
    }
    grpUniq=sort(unique(c(grpUniq1,grpUniq2,grpUniq3,grpUniqSD,grpUniqHTOd)))
    if (length(grpUniq)!=0) {
        colVec=gray(1:length(grpUniq)/length(grpUniq))
        colVec=colVec[!colVec%in%hashInfo$hashColors]
        colVec=colVec[1:length(grpUniq)]
        k=which(grpUniq%in%grpUniq1)
        col_var_info[["snacsRnd1"]]=list(color=c(snacsObj$annHash$hashColors,colVec[k]),level=c(snacsObj$annHash$hashNames,grpUniq[k]))
        k=which(grpUniq%in%grpUniq2)
        col_var_info[["snacsRnd2"]]=list(color=c(snacsObj$annHash$hashColors,colVec[k]),level=c(snacsObj$annHash$hashNames,grpUniq[k]))
        if ("snacsRnd3"%in%names(snacsObj$annCell)) {
            k=which(grpUniq%in%grpUniq3)
            col_var_info[["snacsRnd3"]]=list(color=c(snacsObj$annHash$hashColors,colVec[k]),level=c(snacsObj$annHash$hashNames,grpUniq[k]))
        }
        if ("snacsPlusDoubletD"%in%names(snacsObj$annCell)) {
            k=which(grpUniq%in%grpUniqSD)
            col_var_info[["snacsPlusDoubletD"]]=list(color=c(snacsObj$annHash$hashColors,colVec[k]),level=c(snacsObj$annHash$hashNames,grpUniq[k]))
        }
        if ("HTOdemux"%in%names(snacsObj$annCell)) {
            k=which(grpUniq%in%grpUniqHTOd)
            col_var_info[["hashHTOdemux"]]=list(color=c(snacsObj$annHash$hashColors,colVec[k],"gold","brown"),level=c(snacsObj$annHash$hashNames,grpUniq[k],"Negative","Missing"))
            k=which(col_var_info[["hashHTOdemux"]]$level%in%snacsObj$annCell$HTOdemux)
            col_var_info[["hashHTOdemux"]]=list(color=col_var_info[["hashHTOdemux"]]$color[k],level=col_var_info[["hashHTOdemux"]]$level[k])
        }
        grpUniqP=unique(snacsObj$annCell$patient)
        grpUniqP=grpUniqP[which(!grpUniqP%in%hashInfo$patient)]
        if (length(grpUniqP)!=0) {
            colVec=rainbow(length(grpUniqP)+nrow(hashInfo))
            colVec=gray(1:length(grpUniqP)/length(grpUniqP))
            colVec=colVec[!colVec%in%hashInfo$patColors]
            colVec=colVec[1:length(grpUniqP)]
            col_var_info[["patient"]]=list(color=c(snacsObj$annHash$patColors,colVec),level=c(hashInfo$patient,grpUniqP))
        }
        if ("hashCall_truth"%in%names(snacsObj$annCell)) {
            x=snacsObj$annCell$hashCall_truth; x[is.na(x)]=""
            k=which(patInfo[,1]%in%x)
            col_var_info[["truth"]]=list(color=patInfo[k,2],level=patInfo[k,1])
        }
        for (rId in snacsRndId) {
            colId=paste0("snacsRnd",rId,"_truth_mismatch")
            if (colId%in%names(snacsObj$annCell)) {
                x=snacsObj$annCell[,paste0("snacsRnd",rId,"_truth_mismatch")]; x[is.na(x)]=""
                k=which(patInfo[,1]%in%x)
                col_var_info[[paste0("mismatchRnd",rId)]]=list(color=patInfo[k,2],level=patInfo[k,1])
            }
        }
    } else {
        col_var_info[["snacsRnd1"]]=list(color=c(snacsObj$annHash$hashColors,"cyan2"),level=c(snacsObj$annHash$hashNames,"Doublet"))
        col_var_info[["patient"]]=list(color=c(snacsObj$annHash$patColors,"cyan2"),level=c(gsub("Patient ","pat",snacsObj$annHash$patient),"Doublet"))
    }
    invisible(col_var_info)
}

####################################################################
####################################################################

getTruthCall=function(snacsObj,fName,dirTrueHashCall="../output/accuracy/trueHashCall/") {
    patInfoAll=data.frame(patTruth=c("multiplet","ambiguous",paste0("Sample.",1:nrow(snacsObj$annHash))),hash=c("Multiplet","Ambiguous",snacsObj$annHash$hashNames),pat=c("Multiplet","Ambiguous",sub("Patient ","pat",snacsObj$annHash$patient)))
    annCellInfo=data.frame(colIdOrig=names(snacsObj$annCell),colIdTmp=names(snacsObj$annCell))
    if ("snacsPlusDoubletD"%in%names(snacsObj$annCell)) {
        annCellInfo$colIdTmp[match(c("snacsPlusDoubletD"),names(snacsObj$annCell))]=c("snacsRnd3")
    }
    snacsRndId=1:3; snacsRndId=snacsRndId[which(paste0("snacsRnd",snacsRndId)%in%annCellInfo$colIdTmp)]

    if (file.exists(paste0(dirTrueHashCall,"trueHashCall",fName,".txt"))) {
        hashCall_truth=read.table(paste0(dirTrueHashCall,"trueHashCall",fName,".txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
        snacsObj$annCell$hashCall_truth=hashCall_truth$hashCall_truth[match(snacsObj$annCell$id,hashCall_truth$id)]
        patInfo=patInfoAll[match(sort(unique(snacsObj$annCell$hashCall_truth)),patInfoAll$patTruth),]
        snacsObj$annCell$hashCall_truth=patInfo$pat[match(snacsObj$annCell$hashCall_truth,patInfo$patTruth)]
    } else {
        snacsObj$annCell$hashCall_truth=""
    }

    tbl=matrix(nrow=nrow(snacsObj$annCell),ncol=length(snacsRndId))
    colnames(tbl)=paste0("snacsRnd",snacsRndId,"_truth_mismatch")
    for (rId in snacsRndId) {
        k1=paste0("snacsRnd",rId); k1=annCellInfo$colIdOrig[match(k1,annCellInfo$colIdTmp)]
        colId=paste0("snacsRnd",rId,"_truth_mismatch")
        x=snacsObj$annCell[,k1]
        x[which(!x%in%snacsObj$annHash$hashNames)]="Multiplet"
        x=patInfoAll$pat[match(x,patInfoAll$hash)]
        tbl[,colId]=snacsObj$annCell$hashCall_truth
        tbl[,colId][which(tbl[,colId]==x)]=NA
        tbl[,colId][is.na(snacsObj$annCell$hashCall_truth)]=NA
        #tbl[,colId][which(snacsObj$annCell$hashCall_truth!="Multiplet")]=NA
    }
    snacsObj$annCell=cbind(snacsObj$annCell[,which(!names(snacsObj$annCell)%in%colnames(tbl))],tbl)
    
    invisible(snacsObj)
}

####################################################################
####################################################################
## ------------------------------
## Generate cell annotation tables for accuracy estimation
## Generate heatmaps of multiplet & singlet experiments using only those multiplet SNPs that are found in the singlet experiments
## Parameters
## snacsObj The multi-subject SNACS object for which accuracy is to be calculated.
## snacsRnd SNACS call round on which accuracy will be based on. Default is 1.
## dirData Folder where R objects related to accuracy are to be saved
## dirOutput Folder where output plots and tables are to be saved

#getAnnoForAccuracy <- function(snacsObj,snacsRnd=c(1,2),hashCallMismatch=c(3,2,1),accuracy=F,exptNameSingleSuffix="",dirData="../data/accuracy/",dirTrueHashCall="../output/accuracy/trueHashCall/",dirOutput="../output/accuracy/",writeOutput="ambiguousGeno") {
getAnnoForAccuracy <- function(snacsObj,hashCallMismatch=c(3,2,1),accuracy=F,exptNameSingleSuffix="",dirData="../data/accuracy/",dirTrueHashCall="../output/accuracy/trueHashCall/",dirOutput="../output/accuracy/",writeOutput="ambiguousGeno") {
    #snacsRnd=snacsRnd[1]
    hashCallMismatch=hashCallMismatch[1]
    
    if (F) {
        dirData <- "../data/"
        
        exptNameSingleSuffix="_unfilt"
        fName <- paste0("snacsObj_SNACS5_unfilt.RData")
        fName <- paste0("snacsObj_withDoubletD_SNACS6_unfilt.RData")
        fName <- paste0("snacsObj_withDoubletD_SNACS5_unfilt.RData")
        cat("\n\n----------------- Input file: ",fName,"\n",sep="")
        load(paste0(dirData,fName))
        hashCallMismatch=2; accuracy=F; exptNameSingleSuffix=exptNameSingleSuffix; dirData="../data/accuracy/"; dirTrueHashCall="../output/accuracy/"; dirOutput="../output/accuracy/"; writeOutput="ambiguousGeno"

        dirData <- "../data/"
        numSNP=3
        exptNameSingleSuffix <- paste0("_hashFiltNumSNP",numSNP)
        fName <- paste0("snacsObj_SNACS5",exptNameSingleSuffix,".RData")
        fName <- paste0("snacsObj_withDoubletD_SNACS5",exptNameSingleSuffix,".RData")
        cat("\n\n----------------- Input file: ",fName,"\n",sep="")
        load(paste0(dirData,fName))

        hashCallMismatch=2; accuracy=T; exptNameSingleSuffix=exptNameSingleSuffix; dirData="../data/accuracy/"; dirTrueHashCall="../output/accuracy/trueHashCall/"; dirOutput="../output/accuracy/"; writeOutput="ambiguousGeno"
    }
    
    if (!file.exists(dirData)) dir.create(file.path(dirData))
    if (!file.exists(dirTrueHashCall)) dir.create(file.path(dirTrueHashCall))
    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))
    if (!file.exists(paste0(dirOutput,"heatmapCommonSNP/"))) dir.create(file.path(paste0(dirOutput,"heatmapCommonSNP/")))
    if (!file.exists(paste0(dirOutput,"annCell/"))) dir.create(file.path(paste0(dirOutput,"annCell/")))
    if ("ambiguousGeno"%in%writeOutput) {
        if (!file.exists(paste0(dirOutput,"ambiguousGeno/"))) dir.create(file.path(paste0(dirOutput,"ambiguousGeno/")))
        if (!file.exists(paste0(dirOutput,"snpPosPropMat/"))) dir.create(file.path(paste0(dirOutput,"snpPosPropMat/")))
    }

    dirDataMain="../data/"

    phen=getExptInfoData()
    x=table(phen$run,phen$patient)
    k=apply(x,1,function(x) {sum(x!=0)==1})
    phen=phen[which(phen$run%in%rownames(x)[k]),]
    snacsObj$annHash$exptName=phen$run[match(snacsObj$annHash$patient,phen$patient)]
    
    if ("snacsPlusDoubletD"%in%names(snacsObj$annCell)) {
        names(snacsObj$annCell)[match(c("snacsPlusDoubletD"),names(snacsObj$annCell))]=c("snacsRnd3")
    }
    snacsRndId=1:3; snacsRndId=snacsRndId[which(paste0("snacsRnd",snacsRndId)%in%names(snacsObj$annCell))]
    
    snacsObj$annCell=cbind(snacsObj$annCell,t(round(snacsObj$hashes,3)))
    tbl=t(round(snacsObj$dist2centroidMat,3))
    colnames(tbl)=paste0("dist2centroid_",colnames(tbl))
    snacsObj$annCell=cbind(snacsObj$annCell,tbl)

    patInfoAll=data.frame(patTruth=c("multiplet","ambiguous",paste0("Sample.",1:nrow(snacsObj$annHash))),hash=c("Multiplet","Ambiguous",snacsObj$annHash$hashNames),pat=c("Multiplet","Ambiguous",sub("Patient ","pat",snacsObj$annHash$patient)),exptName=c("Multiplet","Ambiguous",sub("Patient ","pat",snacsObj$annHash$exptName)))

    #fName=paste0("_",snacsObj$exptName,"_",strsplit(snacsObj$exptName,"_")[[1]][1],"snp")
    fName=paste0("_",strsplit(snacsObj$exptName,"_")[[1]][1],exptNameSingleSuffix,"_",strsplit(snacsObj$exptName,"_")[[1]][1],"snp")
    if (file.exists(paste0(dirTrueHashCall,"trueHashCall",fName,".txt"))) {
        accuracy=T
        hashCall_truth=read.table(paste0(dirTrueHashCall,"trueHashCall",fName,".txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T,nrow=-1)
        snacsObj$annCell$hashCall_truth=hashCall_truth$hashCall_truth[match(snacsObj$annCell$id,hashCall_truth$id)]
        patInfo=patInfoAll[match(sort(unique(snacsObj$annCell$hashCall_truth)),patInfoAll$patTruth),]
        snacsObj$annCell$hashCall_truth=patInfo$pat[match(snacsObj$annCell$hashCall_truth,patInfo$patTruth)]
    } else {
        accuracy=F
        snacsObj$annCell$hashCall_truth=""
    }

    tbl=matrix(nrow=nrow(snacsObj$annCell),ncol=length(snacsRndId))
    colnames(tbl)=paste0("snacsRnd",snacsRndId,"_truth_mismatch")
    for (rId in snacsRndId) {
        k1=paste0("snacsRnd",rId); colId=paste0("snacsRnd",rId,"_truth_mismatch")
        x=snacsObj$annCell[,k1]
        x[which(!x%in%snacsObj$annHash$hashNames)]="Multiplet"
        x=patInfoAll$pat[match(x,patInfoAll$hash)]
        tbl[,colId]=snacsObj$annCell$hashCall_truth
        tbl[,colId][which(tbl[,colId]==x)]=NA
        tbl[,colId][is.na(snacsObj$annCell$hashCall_truth) | snacsObj$annCell$hashCall_truth=="Ambiguous"]=NA
        #tbl[,colId][which(snacsObj$annCell$hashCall_truth!="Multiplet")]=NA
    }
    snacsObj$annCell=cbind(snacsObj$annCell,tbl)

    snacsObj0=snacsObj
    exptNameMultInit=paste0(strsplit(snacsObj$exptName,"_")[[1]][1],exptNameSingleSuffix)

    if ("meanQuality"%in%names(snacsObj$annCell)) {
        load(file=paste0(dirDataMain,"snacsObj_init_",exptNameMultInit,".RData"))
        j=match(snacsObj0$annCell$id,snacsObj$annCell$id)
        snacsObj0$annCell$meanQuality=snacsObj$annCell$meanQuality[j]
        snacsObj0$annCell$meanTotalDepth=snacsObj$annCell$meanTotalDepth[j]
    }
    
    snacsObj=snacsObj0

    exptNameMult=snacsObj$exptName
    exptNameSingleVec=paste0(unique(phen$run[phen$patient%in%snacsObj$annHash$patient]),exptNameSingleSuffix)
    exptNameVec=c(exptNameMult,exptNameSingleVec)

    snpName=snacsObj$annSNP$desc
    for (eIdS in 1:length(exptNameSingleVec)) {
        load(file=paste0(dirDataMain,"snacsObj_init_",exptNameSingleVec[eIdS],".RData"))
        snpName=snpName[which(snpName%in%snacsObj$annSNP$desc)]
    }
    if (length(snpName)==0) {
        stop("Accuracy: No commom SNPs among the experiments.")
    }
    
    
    ## ------------------------------
    ## Ambiguous genotypes
    
    if ("ambiguousGeno"%in%writeOutput) {
        snacsObj=snacsObj0
        genotype=as.character(snacsObj$mut[1,])
        if (nrow(snacsObj$mut)>1) {for (i in 2:nrow(snacsObj$mut)) {genotype=paste(genotype,snacsObj$mut[i,])}}
        snacsObj$annCell$genotypeAll=genotype
        id=desc=dirn=samPair=c()
        for (i in 1:nrow(snacsObj$annSNP)) {
            x=strsplit(snacsObj$annSNP$dirn[i],"|",fixed=T)[[1]]
            id=c(id,rep(snacsObj$annSNP$id[i],length(x)))
            desc=c(desc,rep(snacsObj$annSNP$desc[i],length(x)))
            dirn=c(dirn,x)
            samIdThis=sapply(x,function(x) {
                paste(sort(strsplit(gsub("up|down","",x),"_")[[1]]),collapse="_")
            },USE.NAMES=F)
            samPair=c(samPair,samIdThis)
        }
        snpInfo=data.frame(id,dirn,samPair)
        dirnUniq=unique(snpInfo$dirn)
        genoUniq=unique(genotype)
        genoUniqMat=matrix(as.integer(unlist(strsplit(genoUniq," "))),ncol=nrow(snacsObj$annSNP),byrow=T)
        colnames(genoUniqMat)=snacsObj$annSNP$id

        genoUniq=apply(genoUniqMat,1,paste,collapse=" ")
        genoUniqMatCom=genoUniqMat[,which(snacsObj0$annSNP$desc%in%snpName)]
        if (is.matrix(genoUniqMatCom)) {
            genoUniqCom=apply(genoUniqMatCom,1,paste,collapse=" ")
        } else {
            genoUniqMatCom=matrix(genoUniqMatCom,ncol=1)
            colnames(genoUniqMatCom)=colnames(genoUniqMat)[which(snacsObj0$annSNP$desc%in%snpName)]
            genoUniqCom=apply(genoUniqMatCom,1,paste,collapse=" ")
        }
        
        if (F) {
            ii=which(snacsObj$annSNP$desc%in%snpName)
            genotype=as.character(snacsObj$mut[ii[1],])
            if (length(ii)>1) {for (i in 2:length(ii)) {genotype=paste(genotype,snacsObj$mut[ii[i],])}}
            snacsObj$annCell$genotype=genotype
            id=desc=dirn=samPair=c()
            for (i in ii) {
                x=strsplit(snacsObj$annSNP$dirn[i],"|",fixed=T)[[1]]
                id=c(id,rep(snacsObj$annSNP$id[i],length(x)))
                desc=c(desc,rep(snacsObj$annSNP$desc[i],length(x)))
                dirn=c(dirn,x)
                samIdThis=sapply(x,function(x) {
                    paste(sort(strsplit(gsub("up|down","",x),"_")[[1]]),collapse="_")
                },USE.NAMES=F)
                samPair=c(samPair,samIdThis)
            }
            snpInfoCom=data.frame(id,dirn,samPair)
            #dirnAmb=dirnUniq[!dirnUniq%in%snpInfoCom$dirn]
        }
        
        dirnAmb=snpInfo$samPair

        if (length(dirnAmb)!=0) {
            if (F) {
                genoAmpSamInfo=data.frame(exptName=rep(strsplit(snacsObj$exptName,"_")[[1]][1],length(dirnAmb)),dirn=dirnAmb,samId1=dirnAmb,samId2=dirnAmb)
                for (k in 1:length(dirnAmb)) {
                    #x=strsplit(gsub("up|down","",dirnAmb[k]),"_")[[1]]
                    x=strsplit(dirnAmb[k],"_")[[1]]
                    samIdThis=patInfoAll$exptName[match(x,patInfoAll$hash)]
                    genoAmpSamInfo$samId1[k]=samIdThis[1]
                    genoAmpSamInfo$samId2[k]=samIdThis[2]
                }
            }
            genoInfo=data.frame(exptName=rep(strsplit(snacsObj$exptName,"_")[[1]][1],length(genoUniqCom)),genotype=genoUniqCom)
        } else {
            #genoAmpSamInfo=genoInfo=NULL
            genoInfo=NULL
        }
    }


    ## ------------------------------
    #snpMat=matrix(nrow=length(snpName),ncol=length(exptNameVec),dimnames=list(snacsObj$annSNP$id[match(snpName,snacsObj$annSNP$desc)],sapply(exptNameVec,function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)))
    snpMat=matrix(nrow=length(snpName),ncol=length(exptNameVec))
    colnames(snpMat)=exptNameVec
    
    ## Create annCell tables containing observed "genotype" column for single and multiple patient experiments
    for (eId in 1:length(exptNameVec)) {
        if (exptNameVec[eId]==exptNameMult) {
            ## Create SNACS object for multi-patient experient with common SNPs
            snacsObj=snacsObj0
            
            if (T) {
                
                ## Get unimputed mutation
                snacsObj$mutUnimp=snacsObj$mut
                i=match(snacsObj$annSNP$id,rownames(snacsObj$missing))
                j=match(snacsObj$annCell$id,colnames(snacsObj$missing))
                snacsObj$missing=snacsObj$missing[i,j]
                snacsObj$mutUnimp[snacsObj$missing]=9
                if (length(i)==1) {
                    snacsObj$missing=matrix(snacsObj$missing,nrow=1)
                }

                ii=match(snpName,snacsObj$annSNP$desc)
                ## Get genotype of SNPs not in single experiments for each cell in the experiment
                if (1%in%ii) genotype=rep("-",ncol(snacsObj$mut)) else genotype=as.character(snacsObj$mut[1,])
                if (nrow(snacsObj$mut)>1) {
                    for (i in 2:nrow(snacsObj$mut)) {
                        if (i%in%ii) genotype=paste(genotype,rep("-",ncol(snacsObj$mut))) else genotype=paste(genotype,snacsObj$mut[i,])
                    }
                }
                snacsObj$annCell$genotypeNA=genotype
                ## Get genotype of SNPs not in single experiments for each cell in the experiment for unimputed data
                if (1%in%ii) genotype=rep("-",ncol(snacsObj$mutUnimp)) else genotype=as.character(snacsObj$mutUnimp[1,])
                if (nrow(snacsObj$mutUnimp)>1) {
                    for (i in 2:nrow(snacsObj$mutUnimp)) {
                        if (i%in%ii) genotype=paste(genotype,rep("-",ncol(snacsObj$mutUnimp))) else genotype=paste(genotype,snacsObj$mutUnimp[i,])
                    }
                }
                snacsObj$annCell$genoUnimpNA=genotype
            }

            i=match(snpName,snacsObj$annSNP$desc)
            #j=match(annCellAll$id,snacsObj$annCell$id)
            j=1:nrow(snacsObj$annCell)
            snacsObj$mut=snacsObj$mut[i,j]
            snacsObj$annSNP=snacsObj$annSNP[i,]
            snacsObj$annCell=snacsObj$annCell[j,]
            snacsObj$hashes=snacsObj$hashes[,j]
            if (length(i)==1) {
                snacsObj$mut=matrix(snacsObj$mut,nrow=1)
            }
            
            ## Get unimputed mutation
            if (F) {
                snacsObj$mutUnimp=snacsObj$mut
                i=match(snacsObj$annSNP$id,rownames(snacsObj$missing))
                j=match(snacsObj$annCell$id,colnames(snacsObj$missing))
                snacsObj$missing=snacsObj$missing[i,j]
                snacsObj$mutUnimp[snacsObj$missing]=9
                if (length(i)==1) {
                    snacsObj$missing=matrix(snacsObj$missing,nrow=1)
                }
            } else {
                snacsObj$mutUnimp=snacsObj$mutUnimp[i,j]
                if (length(i)==1) {
                    snacsObj$mutUnimp=matrix(snacsObj$mutUnimp,nrow=1)
                }
            }

            snacsObj <- getHTOdemuxCall(snacsObj=snacsObj)

            hashUniq=unique(c(snacsObj$annHash$hashNames,sort(unique(snacsObj$annCell$snacsRnd1))))
            patUniq=gsub("Patient ","pat",snacsObj$annHash$patient)
            snacsObj$annCell$patient=patUniq[match(snacsObj$annCell$snacsRnd1,hashUniq)]

            snacsObj$annCell$patient=ifelse(sum(!duplicated(snacsObj$annHash$patient))==1,snacsObj$annHash$patient[1],"Multiplet")
            j=match(snacsObj$annCell$snacsRnd1,snacsObj$annHash$hashNames); j1=which(!is.na(j)); j2=j[j1]
            snacsObj$annCell$patient[j1]=snacsObj$annHash$patient[j2]
            snacsObj$annCell$patient=gsub("Patient ","pat",snacsObj$annCell$patient)
            
            tbl=matrix(nrow=nrow(snacsObj$annCell),ncol=length(snacsRndId))
            colnames(tbl)=paste0("patientRnd",snacsRndId)
            for (rId in snacsRndId) {
                k1=paste0("snacsRnd",rId); colId=paste0("patientRnd",rId)
                tbl[,colId]=ifelse(sum(!duplicated(snacsObj$annHash$patient))==1,snacsObj$annHash$patient[1],"Multiplet")
                j=match(snacsObj$annCell[,k1],snacsObj$annHash$hashNames); j1=which(!is.na(j)); j2=j[j1]
                tbl[j1,colId]=snacsObj$annHash$patient[j2]
                tbl[,colId]=gsub("Patient ","pat",tbl[,colId])
            }
            snacsObj$annCell=cbind(snacsObj$annCell,tbl)
            
            if ("ambiguousGeno"%in%writeOutput) {
                snpMat[match(snpName,snacsObj$annSNP$desc),which(colnames(snpMat)==snacsObj$exptName)]=apply(snacsObj$mut,1,function(x) {mean(x==1,na.rm=T)})
            }

            col_var_info=getColVarInfo(snacsObj,snacsRndId=snacsRndId)
        } else {
            ## Create SNACS object for single-patient experient with common SNPs
            snpName2=snpName
            load(file=paste0(dirDataMain,"snacsObj_imp_",exptNameVec[eId],".RData"))
            snpName2=snacsObj$annSNP$desc

            load(file=paste0(dirDataMain,"snacsObj_init_",exptNameVec[eId],".RData"))
            i=which(snacsObj$annSNP$desc%in%c(snpName2,snacsObj0$annSNP$desc))
            snacsObj$mut=snacsObj$mut[i,]
            snacsObj$annSNP=snacsObj$annSNP[i,]
            #table(snpName%in%snacsObj$annSNP$desc)
            
            #snacsObj=imputeMissingMutations(snacsObj,proportionMissingPerCell=1,proportionMissingPerSNP=1,proportionMutatedPerCell=c(0,1),proportionMutatedPerSNP=c(0,1))
            snacsObj=filterData(snacsObj,proportionMissingPerCell=1,proportionMissingPerSNP=1,proportionMutatedPerCell=c(0,1),proportionMutatedPerSNP=c(0,1))
            snacsObj=imputeMissingMutations(snacsObj)

            i=match(snpName,snacsObj$annSNP$desc)
            snacsObj$mut=snacsObj$mut[i,]
            snacsObj$annSNP=snacsObj$annSNP[i,]
            if (length(i)==1) {
                snacsObj$mut=matrix(snacsObj$mut,nrow=1)
            }

            ## Get unimputed mutation
            snacsObj$mutUnimp=snacsObj$mut
            i=match(snacsObj$annSNP$id,rownames(snacsObj$missing))
            j=match(snacsObj$annCell$id,colnames(snacsObj$missing))
            snacsObj$missing=snacsObj$missing[i,j]
            snacsObj$mutUnimp[snacsObj$missing]=9
            if (length(i)==1) {
                snacsObj$missing=matrix(snacsObj$missing,nrow=1)
            }

            snacsObj$annCell=cbind(snacsObj$annCell,t(round(snacsObj$hashes,3)))

            snacsObj$annCell$patient=snacsObj$annHash$patient[1]
            snacsObj$annCell$patient=gsub("Patient ","pat",snacsObj$annCell$patient)
            snacsObj$annCell$snacsRnd1="Multiplet"
            
            if ("ambiguousGeno"%in%writeOutput) {
                if (length(dirnAmb)!=0) {
                    genotype=as.character(snacsObj$mut[1,])
                    if (nrow(snacsObj$mut)>1) {for (i in 2:nrow(snacsObj$mut)) {genotype=paste(genotype,snacsObj$mut[i,])}}
                    genoUniq=unique(genotype)
                    genoInfo=rbind(genoInfo,data.frame(exptName=rep(strsplit(snacsObj$exptName,"_")[[1]][1],length(genoUniq)),genotype=genoUniq))
                }
                snpMat[match(snpName,snacsObj$annSNP$desc),which(colnames(snpMat)==snacsObj$exptName)]=apply(snacsObj$mut,1,function(x) {mean(x==1,na.rm=T)})
            }
            
            col_var_info=list()
            col_var_info[["patient"]]=list(color=c(snacsObj$annHash$patColors[1]),level=c(gsub("Patient ","pat",snacsObj$annHash$patient[1])))
        }
        
        col_var_info[["meanQuality"]]=list(limit=c(10,90))
        col_var_info[["meanTotalDepth"]]=list(limit=c(20,120))

        ## Get genotype for each cell in the experiment
        genotype=as.character(snacsObj$mut[1,])
        if (nrow(snacsObj$mut)>1) {for (i in 2:nrow(snacsObj$mut)) {genotype=paste(genotype,snacsObj$mut[i,])}}
        snacsObj$annCell$genotype=genotype

        ## Get genotype for each cell in the experiment for unimputed data
        genotype=as.character(snacsObj$mutUnimp[1,])
        if (nrow(snacsObj$mutUnimp)>1) {for (i in 2:nrow(snacsObj$mutUnimp)) {genotype=paste(genotype,snacsObj$mutUnimp[i,])}}
        snacsObj$annCell$genoUnimputed=genotype

        ## ------------------------------
        ## Create heatmap
        annCell=as.data.frame(t(snacsObj$hashes[nrow(snacsObj$hashes):1,]))
        nm=rev(snacsObj$annHash$hashNames)
        if (exptNameVec[eId]==exptNameMult) {
            #annCell=snacsObj$annCell[,c("genotype","patient","snacsRnd1","clustSplitRnd1","snacsRnd2","clustSplitRnd2","HTOdemux")]
            #names(annCell)=c("genotype","patient","hashRnd1","clustSplitRnd1","hashRnd2","clustSplitRnd2","hashHTOdemux")
            #annCell=cbind(snacsObj$annCell[,c("meanQuality","meanTotalDepth","snacsRnd1","clustSplitRnd1",paste0("snacsRnd",hashCallMismatch,"_truth_mismatch"),"hashCall_truth","snacsRnd2","clustSplitRnd2")],annCell)
            #nm=c(c("quality","depth","hashRnd1","clustSplitRnd1","mismatch","truth","hashRnd2","clustSplitRnd2"),nm)
            #annCell=cbind(snacsObj$annCell[,c("meanQuality","meanTotalDepth","snacsRnd1","clustSplitRnd1","snacsRnd3",paste0("snacsRnd",hashCallMismatch,"_truth_mismatch"),"hashCall_truth","snacsRnd2","clustSplitRnd2")],annCell)
            #nm=c(c("quality","depth","hashRnd1","clustSplitRnd1","hashRnd3","mismatch","truth","hashRnd2","clustSplitRnd2"),nm)

            for (rId in snacsRndId) {
                colId=paste0(c("snacsRnd","clustSplitRnd"),rId); colId=c(colId,paste0("snacsRnd",rId,"_truth_mismatch"))
                colName=paste0(c("hashRnd","clustSplitRnd","mismatchRnd"),rId)
                k=which(colId%in%names(snacsObj$annCell))
                annCell=cbind(snacsObj$annCell[,colId[k]],annCell)
                nm=c(colName[k],nm)
            }
            if ("doubletD"%in%names(snacsObj$annCell)) {
                annCell=cbind(doubletD=snacsObj$annCell$doubletD,annCell)
                nm=c("doubletD",nm)
            }
            annCell=cbind(hashCall_truth=snacsObj$annCell$hashCall_truth,annCell)
            nm=c("truth",nm)
            #for (k in grep("truth",names(annCell))) annCell[which(is.na(annCell[,k]) | annCell[,k]=="Ambiguous"),k]=""
            for (k in grep("truth",names(annCell))) annCell[is.na(annCell[,k]),k]=""

            if (F) {
                #patInfo=cbind(sort(unique(annCell$hashCall_truth)),c("Multiplet",sub("Patient ","pat",snacsObj0$annHash$patient)))
                #j=match(annCell$hashCall_truth,patInfo[,1])
                #annCell$hashCall_truth=patInfo[j,2]
                annCell$hashCall_truth[is.na(annCell$hashCall_truth)]=""
                k=paste0("snacsRnd",hashCallMismatch,"_truth_mismatch")
                annCell[is.na(annCell[,k]),k]=""
            }
        } else {
            #annCell=snacsObj$annCell[,c("genotype","patient")]
            #names(annCell)=c("genotype","patient")
            colId=c("meanQuality","meanTotalDepth","patient")
            colName=c("quality","depth","patient")
            k=which(colId%in%names(snacsObj$annCell))
            if (length(k)!=0) {
                annCell=cbind(snacsObj$annCell[,colId[k]],annCell)
                nm=c(colName[k],nm)
            }
        }
        colNameInfo=data.frame(colId=names(annCell),colName=nm)
        colNameInfo$colName=sub("hashRnd3","snacsPlusDoubletD",colNameInfo$colName)
        names(annCell)=colNameInfo$colName
        k=match(names(col_var_info),colNameInfo$colId); k1=which(!is.na(k)); k2=k[k1]
        names(col_var_info)[k1]=colNameInfo$colName[k2]
        snacsObj$col_var_info=col_var_info
        
        fName=paste0("_",snacsObj$exptName,"_",strsplit(exptNameMult,"_")[[1]][1],"snp")
        #save(snacsObj,file=paste0(dirData,"snacsObj",fName,".RData"))

        header=paste0(snacsObj$exptName," (",snacsObj0$exptName," SNPs)")
        nameRow=paste0("SNP-",nrow(snacsObj$annSNP):1)
        if (exptNameVec[eId]==exptNameMult) {
            #nameRow=paste0(nameRow," ",rev(snacsObj$annSNP$desc))
            nameRow=paste0(nameRow," ",snacsObj$annSNP$desc)
        }
        if (nrow(snacsObj$annSNP)==1) {
            x=rbind(snacsObj$mut,snacsObj$mut)
            nameRow=c(nameRow,"")
        } else {
            x=snacsObj$mut
        }
        plot_info=list(cexRow=.8)

        col_dend=T; row_dend=F
        zlm=c(0,1); heatmap_color=c("navy","gray85","orangered1")
        input_legend=F
        #png(paste0("heatmap",fName,".png"))
        pdf(paste0(dirOutput,"heatmapCommonSNP/heatmap",fName,".pdf"))
        clustObj=heatmap4::generate_heatmap(x=x,col_clust=snacsObj$hclustObj_bestSNPs,col_dend=T,row_dend=row_dend,col_info=annCell,col_anno=T,zlm=zlm,heatmap_color=heatmap_color,col_var_info=snacsObj$col_var_info,h_title=header,row_lab=T,row_lab_vtr=nameRow,input_legend=input_legend,plot_info=plot_info,densColor=30)
        dev.off()

        ## ------------------------------
        ## Write out cell-level table having genotype information
        #colId=c("id","patient","meanQuality","meanTotalDepth","snacsRnd1","clustSplitRnd1","snacsRnd3","hashCall_truth",paste0("snacsRnd",snacsRndId,"_truth_mismatch"),"snacsRnd2","clustSplitRnd2","genotype","genoUnimputed")
        tbl=snacsObj$annCell[clustObj$colClust$order,]
        tbl=tbl[,which(!names(tbl)%in%"hashRnd1")]
        if (all(tbl$snacsRnd1=="Multiplet")) tbl=tbl[,which(!names(tbl)%in%"snacsRnd1")]
        k=match(c("snacsRnd3","snacsRnd3_truth_mismatch"),names(tbl)); k1=which(!is.na(k))
        if (length(k1)!=0) names(tbl)[k[k1]]=c("snacsPlusDoubletD","snacsPlusDoubletD_truth_mismatch")[k1]
        write.table(tbl,file=paste0(dirOutput,"annCell/annCell",fName,".txt"),col.names=T,row.names=F, sep="\t",quote=F)
        
        if (accuracy & exptNameVec[eId]==exptNameMult) {
            if (hashCallMismatch==3) {colId="snacsRnd3"
            } else if (hashCallMismatch==2) {colId="snacsRnd2"
            } else {colId="snacsRnd1"}
            patInfo=patInfoAll[match(c(snacsObj$annHash$hashNames,"Multiplet","Ambiguous"),patInfoAll$hash),]
            #patInfo=cbind(c(snacsObj$annHash$hashNames,"Multiplet"),c(sub("Patient ","pat",snacsObj$annHash$patient),"Multiplet"))
            x0=snacsObj$annCell$hashCall_truth
            x0[which(x0=="Ambiguous")]=NA
            x=snacsObj$annCell[,colId]
            x[which(!x%in%snacsObj$annHash$hashNames)]="Multiplet"
            x=patInfo$pat[match(x,patInfo$hash)]
            y=data.frame(truth=snacsObj$annCell$hashCall_truth,snacs=x)
            y=y[which(!y$truth%in%c("Multiplet","Ambiguous")),]
            #cat("\nAccuracy based on round ",hashCallMismatch,": ",round(mean(snacsObj$annCell$hashCall_truth==x,na.rm=T),2),"\n\n",sep="")
            cat("\nSens / spec / acc based on round ",hashCallMismatch,": ",
                round(sum(snacsObj$annCell$hashCall_truth=="Multiplet" & x=="Multiplet",na.rm=T)/sum(snacsObj$annCell$hashCall_truth=="Multiplet",na.rm=T),2)," / ",
                round(mean(y$truth==y$snacs,na.rm=T),2)," / ",
                round(mean(x0==x,na.rm=T),2),"\n\n",sep="")
            #print(table(truth=(snacsObj$annCell$hashCall_truth!="Multiplet"),snacs=(x!="Multiplet"),exclude=NULL))
            print(table(truth=snacsObj$annCell$hashCall_truth,snacs=x,exclude=NULL))
        }
    }
    
    if ("ambiguousGeno"%in%writeOutput) {
        if (length(dirnAmb)!=0) {
            snacsObj=snacsObj0
            
            k=which(genoInfo$exptName==strsplit(exptNameMult,"_")[[1]][1])
            tmpC=rep("",length(k))
            genoAmpInfo=cbind(genoInfo[k,],samId1=tmpC,samId2=tmpC)
            grpUniq=unique(genoInfo$exptName[-k])
            for (gId1 in 1:(length(grpUniq)-1)) {
                j1=which(genoInfo$exptName==grpUniq[gId1])
                for (gId2 in (gId1+1):length(grpUniq)) {
                    j2=which(genoInfo$exptName==grpUniq[gId2])
                    i=which(genoAmpInfo$genotype%in%intersect(genoInfo$genotype[j1],genoInfo$genotype[j2]))
                    if (length(i)!=0) {
                        genoAmpInfo$samId1[i]=grpUniq[gId1]
                        genoAmpInfo$samId2[i]=grpUniq[gId2]
                    }
                }
            }
            if (any(genoAmpInfo$samId1!="")) {
                genoAmpInfo=genoAmpInfo[!duplicated(paste(genoAmpInfo$genotype,genoAmpInfo$samId1,genoAmpInfo$samId2)) & genoAmpInfo$samId1!="",]
                #genoAmpInfo=genoAmpInfo[genoAmpInfo$samId1!="",]
                fName=paste0("_",snacsObj$exptName,"_",strsplit(exptNameMult,"_")[[1]][1],"snp")
                write.table(genoAmpInfo,file=paste0(dirOutput,"ambiguousGeno/ambiguousGeno",fName,".txt"),col.names=T,row.names=F, sep="\t",quote=F)
            }
        }
        colnames(snpMat)=sapply(colnames(snpMat),function(x) {strsplit(x,"_")[[1]][1]},USE.NAMES=F)
        fName=paste0("_",snacsObj$exptName,"_",strsplit(exptNameMult,"_")[[1]][1],"snp")
        tbl=cbind(id=snacsObj0$annSNP$id[match(snpName,snacsObj0$annSNP$desc)],desc=snpName,as.data.frame(round(snpMat,5))); tbl=tbl[nrow(tbl):1,]
        write.table(tbl,file=paste0(dirOutput,"snpPosPropMat/snpPosPropMat",fName,".txt"),col.names=T,row.names=F, sep="\t",quote=F)
    }
}

## ------------------------------
clusterSampleWithAntibodyData.my=function(snacsObj,backgndThreshold=0.95,cellProportionAboveBackgnd=0.5,cellProportionBelowBackgndMode=0.6,cellProportionForModeDetection=0.75,hashThreshold=0.5) {
    ## --------------------------------------
    dirOutput="../output/hashPairPlot/"
    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))

    ## --------------------------------------
    ## --------------------------------------
    #backgndThreshold=0.95; cellProportionBelowBackgndMode=0.6; cellProportionForModeDetection=0.75; cellProportionAboveBackgnd=0.5
    backgndThresRnd2=0.75
    
    #minClustSize=2; clustComparePValue=10^-5; maxClustSampleSize=Inf; clustCompareMethod=c("t","hotelling")
    #clustCompareMethod[clustCompareMethod[1]]
    
    hashMat=t(snacsObj$hashes)
    ## --------------------------------------
    ## --------------------------------------

    hashBackgnd=matrix(nrow=4,ncol=nrow(snacsObj$annHash),dimnames=list(c("mean","sd","thres","thresRnd2"),snacsObj$annHash$hashNames))
    for (k in 1:ncol(hashMat)) {
        x=stats::density(hashMat[,k],bw="SJ",na.rm=T)
        xlim=range(x$x); ylim=range(x$y); ylim=NULL; ylim=c(0,1)
        xlim=max(abs(range(x$x))); xlim=c(-xlim,xlim)
        k1=which.max(x$y)
        if (stats::quantile(hashMat[,k],probs=cellProportionBelowBackgndMode)<x$x[k1]) k1=which.max(x$y[1:stats::quantile(1:(k1-1),probs=cellProportionForModeDetection)])
        xx=x$x[k1]
        x2=x$x[which(x$x<xx)]-xx; x2=c(x2,-x2); x2=x2+xx
        hashBackgndData=x2
        hashBackgndECDF=stats::ecdf(hashBackgndData)
        x2=hashMat[which(hashMat[,k]<xx),k]-xx; x2=c(x2,-x2); x2=x2+xx
        ecdfbg=1-hashBackgndECDF(x2)
        hashBackgnd[,k]=c(mean(x2),stats::sd(x2),stats::quantile(x2,probs=c(backgndThreshold,backgndThresRnd2)))
    }
    
    hashCallThis=rep("",nrow(snacsObj$annCell))
    hashMatThis=hashMat
    j=1:length(hashCallThis)
    for (sampleId in snacsObj$annHash$hashNames) {
        if (backgndThreshold==0) {
            jj=which(hashMatThis[,sampleId]>=cellProportionAboveBackgnd)
        } else {
            jj=which(hashMatThis[,sampleId]>hashBackgnd["thres",sampleId])
        }
        if (length(jj)!=0) hashCallThis[jj]=paste0(hashCallThis[jj],"_",sampleId)
    }
    hashCallThis=sub("^_", "",hashCallThis)
    
    snacsObj$annCell$hashClust=hashCallThis
    
    if (F) {
        fName=paste0(backgndThreshold,"_",cellProportionAboveBackgnd,"_",cellProportionBelowBackgndMode,"_",cellProportionForModeDetection)
        header=paste0(snacsObj$exptName,"\nBgnd: ",backgndThreshold," thres, ",cellProportionAboveBackgnd," above, ",cellProportionBelowBackgndMode," below, ",cellProportionForModeDetection," mode")
        header=paste0(snacsObj$exptName,"\nBgnd: ",backgndThreshold," thres, ",cellProportionAboveBackgnd," > prop cell, ",cellProportionBelowBackgndMode," below, ",cellProportionForModeDetection," mode")
        #grDevices::pdf(paste0(dirOutput,"scatterPlot_hashPair_",snacsObj$exptName,"_",fName,".pdf"))
        for (samId1 in 1:(nrow(snacsObj$annHash)-1)) {
            j1=which(snacsObj$annCell$hashClust==snacsObj$annHash$hashNames[samId1])
            for (samId2 in (samId1+1):nrow(snacsObj$annHash)) {
                j2=which(snacsObj$annCell$hashClust==snacsObj$annHash$hashNames[samId2])
                graphics::plot(snacsObj$hashes[samId1,],snacsObj$hashes[samId2,],main=header,xlab=paste0(snacsObj$annHash$hashNames[samId1],": Hash"),ylab=paste0(snacsObj$annHash$hashNames[samId2],": Hash"),pch=16)
                graphics::points(snacsObj$hashes[samId1,j1],snacsObj$hashes[samId2,j1],pch=16,col=snacsObj$annHash$hashColors[samId1])
                graphics::points(snacsObj$hashes[samId1,j2],snacsObj$hashes[samId2,j2],pch=16,col=snacsObj$annHash$hashColors[samId2])
            }
        }
        #grDevices::dev.off()
    }

    invisible(snacsObj)
}
####################################################################
####################################################################
