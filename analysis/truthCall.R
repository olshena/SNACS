## SNACS_Manuscript

####################################################################
####################################################################
## Get unique logical genotype-combinations from all genotype-combinations from an experiment
getUniqueGenoComb=function(exptName,thresGenoCall=0,thresGenoNA=NULL,dirData="../output/accuracy/snpPosPropMat/") {
    #exptName="SNACS15"; thresGenoCall=thresGenoCall; dirData="../output/accuracy/snpPosPropMat/"
    
    cat("\n------------ getUniqueGenoComb: ",exptName," --------------\n",sep="")
    timeStamp=Sys.time()
    
    fName=paste0(dirData,"snpPosPropMat_",exptName,exptNameSuffix,"_",exptName,"snp.txt")
    if (file.exists(fName)) {
        tbl=read.delim(fName,sep="\t",header=TRUE,as.is=TRUE)
        snpMat=as.matrix(tbl[nrow(tbl):1,4:ncol(tbl)])
        snpMat[snpMat<=thresGenoCall]=0
        snpMat[snpMat>=(1-thresGenoCall)]=1
        snpMat[snpMat>thresGenoCall & snpMat<(1-thresGenoCall)]=0.5
        if (!is.null(thresGenoNA)) {
            snpMat[snpMat>thresGenoCall & snpMat<thresGenoNA]=NA
            snpMat[snpMat<(1-thresGenoCall) & snpMat>(1-thresGenoNA)]=NA
        }
        genoUniq=list()
        for (k in 1:ncol(snpMat)) {
            x=snpMat[!is.na(snpMat[,k]),k]
            if (length(x)==0) {
                y=NULL
            } else {
                i=which(x==0.5)
                if (length(i)==0) {
                    y=paste0(x,collapse="")
                } else {
                    y=c()
                    n = length(i)
                    z <- rep(0, n)
                    z=do.call(rbind, lapply(0:n, function(i) t(apply(combn(1:n,i), 2, function(k) {z[k]=1;z}))))
                    x1=rep(x,nrow(z))
                    for (p in 1:nrow(z)) {
                        x1=x
                        x1[i]=z[p,]
                        y=c(y,paste0(x1,collapse=""))
                    }

                }
            }
            genoUniq[[k]]=unique(y)
        }
        names(genoUniq)=colnames(snpMat)
    } else {
        genoUniq=NULL
    }
    
    timeStamp=c(timeStamp,Sys.time())
    print(timeStamp)
    print(diff(timeStamp))
    
    invisible(genoUniq)
}

getAmbiguousGenotype=function(exptName,thresGenoCall=0,thresGenoNA=NULL,dirData="../output/accuracy/snpPosPropMat/") {
    cat("\n------------ getAmbiguousGenotype: ",exptName," --------------\n",sep="")
    timeStamp=Sys.time()
    
    fName=paste0(dirData,"snpPosPropMat_",exptName,exptNameSuffix,"_",exptName,"snp.txt")
    if (file.exists(fName)) {
        tbl=read.delim(fName,sep="\t",header=TRUE,as.is=TRUE)
        snpMat=as.matrix(tbl[nrow(tbl):1,4:ncol(tbl)])
        snpMat[snpMat<=thresGenoCall]=0
        snpMat[snpMat>=(1-thresGenoCall)]=1
        snpMat[snpMat>thresGenoCall & snpMat<(1-thresGenoCall)]=0.5
        if (!is.null(thresGenoNA)) {
            snpMat[snpMat>thresGenoCall & snpMat<thresGenoNA]=NA
            snpMat[snpMat<(1-thresGenoCall) & snpMat>(1-thresGenoNA)]=NA
        }
        genoUniq=list()
        for (k in 1:ncol(snpMat)) {
            x=snpMat[!is.na(snpMat[,k]),k]
            if (length(x)==0) {
                y=NULL
            } else {
                i=which(x==0.5)
                if (length(i)==0) {
                    y=paste0(x,collapse="")
                } else {
                    y=c()
                    n = length(i)
                    z <- rep(0, n)
                    z=do.call(rbind, lapply(0:n, function(i) t(apply(combn(1:n,i), 2, function(k) {z[k]=1;z}))))
                    x1=rep(x,nrow(z))
                    for (p in 1:nrow(z)) {
                        x1=x
                        x1[i]=z[p,]
                        y=c(y,paste0(x1,collapse=""))
                    }

                }
            }
            genoUniq[[k]]=unique(y)
        }
        names(genoUniq)=colnames(snpMat)
        amb=NULL
        for (k1 in 1:(length(genoUniq)-1)) {
            for (k2 in (k1+1):length(genoUniq)) {
                x=intersect(genoUniq[[k1]],genoUniq[[k2]])
                if (length(x)!=0) {
                    amb=rbind(amb,data.frame(samId1=names(genoUniq)[k1],samId2=names(genoUniq)[k2],genotype=x))
                }
            }
        }
        if (!is.null(amb)) {
            amb=amb[!duplicated(paste(amb$samId1,amb$samId2,amb$genotype)),]
            exptNameS=names(genoUniq)
            samNames=paste0("Sample.",as.integer(as.factor(exptNameS)))
            amb$samId1=samNames[match(amb$samId1,exptNameS)]
            amb$samId2=samNames[match(amb$samId2,exptNameS)]
        }
    } else {
        amb=NULL
    }
    
    timeStamp=c(timeStamp,Sys.time())
    print(timeStamp)
    print(diff(timeStamp))

    invisible(amb)
}


makeHashPairScatterPlot=function(exptName,exptNameSuffix,dirData=paste0(dirOutput,"/hashPairPlot/")) {
    load(paste0("../data/snacsObj_",exptName,exptNameSuffix,".RData"))
    x=NULL
    if (exptName=="SNACS5") {x=cbind(dat5,foo=foo5);
    } else if (exptName=="SNACS6") {x=cbind(dat6,foo=foo6)
    } else if (exptName=="SNACS7") {x=cbind(dat7,foo=foo7)}
    snacsObj$annCell$hashCall_truth[match(rownames(x),snacsObj$annCell$id)]=x$foo
    j=match(snacsObj$annCell$hashCall_truth,c(paste0("Sample.",1:nrow(snacsObj$annHash)),"multiplet"))
    j1=which(!is.na(j)); j2=j[j1]
    snacsObj$annCell$hashCall_truth=""
    snacsObj$annCell$hashCall_truth[j1]=c(snacsObj$annHash$hashNames,"Multiplet")[j2]
    colVec=rep("gray",nrow(snacsObj$annCell)); colVec[j1]=c(snacsObj$annHash$hashColor,"black")[j2]
    fName=paste0("")
    header=paste0(snacsObj$exptName,": Truth calls")
    #dirData=paste0("../output/accuracy/truthCall",ifelse(truthCallType=="manual","Manual",paste0("Auto_thres",thresGenoCall)),"/hashPairPlot/")
    if (!file.exists(dirData)) dir.create(file.path(dirData))
    pdf(paste0(dirData,"scatterPlot_hashPair_truthCall_",snacsObj$exptName,".pdf"))
    #j=which(snacsObj$annCell$hashCall_truth=="TS.2" & snacsObj$annCell$snacsRnd2=="TS.4")
    for (samId1 in 1:(nrow(snacsObj$annHash)-1)) {
        #j1=which(snacsObj$annCell$hashClust==snacsObj$annHash$hashNames[samId1])
        for (samId2 in (samId1+1):nrow(snacsObj$annHash)) {
            #j2=which(snacsObj$annCell$hashClust==snacsObj$annHash$hashNames[samId2])
            plot(snacsObj$hashes[samId1,],snacsObj$hashes[samId2,],main=header,xlab=paste0(snacsObj$annHash$hashNames[samId1],": Hash"),ylab=paste0(snacsObj$annHash$hashNames[samId2],": Hash"),pch=16,col=colVec)
            #points(snacsObj$hashes[samId1,j1],snacsObj$hashes[samId2,j1],pch=16,col=snacsObj$annHash$hashColors[samId1])
            #points(snacsObj$hashes[samId1,j2],snacsObj$hashes[samId2,j2],pch=16,col=snacsObj$annHash$hashColors[samId2])
            #points(snacsObj$hashes[samId1,j],snacsObj$hashes[samId2,j],pch="x",col="green3")
        }
    }
    dev.off()
}

####################################################################
####################################################################

makeTruthCall=function(exptNameSuffix,thresGenoCall,thresGenoNA=NULL,dirData="../output/accuracy/annCell/",dirOutput="../output/accuracy/trueHashCall/",saveFlag=F,loadFlag=T,maxSam=4) {
    #maxSam=8; saveFlag=T; loadFlag=T; exptNameSuffix=exptNameSuffix; thresGenoCall=thresGenoCall; thresGenoNA=thresGenoNA; dirData="../output/accuracy/annCell/"; dirOutput="../output/accuracy/trueHashCall/"
    
    
    #truthCallType="manual" ## NOT used. For hand curated singlet truth calls
    truthCallType="automatic"
    
    timeStamp=Sys.time()
    print(timeStamp)

    cat("\n\n---------------------------------\n")
    cat("------------------",exptNameSuffix," thres ",thresGenoCall," ---------------\n",sep="")
    cat("---------------------------------\n\n")

    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))

    
    #dirOutput=paste0("../output/accuracy/truthCall",ifelse(truthCallType=="manual","Manual",paste0("Auto_thres",thresGenoCall)),"/trueHashCall/")

    if (F) {
    ####################################################################
    ## SNACS5

    exptMultId=5
    exptSingleId=c(1,2)
    exptIdVec=c(exptMultId,exptSingleId)

    nVec=rep(NA,length(exptIdVec))
    for (eId in 1:length(exptIdVec)) {
        tbl=read.table(paste0(dirData,"annCell_SNACS",exptIdVec[eId],exptNameSuffix,"_SNACS",exptMultId,"snp.txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        for (k in c("genotype","genoUnimputed")) tbl[,k]=as.character(tbl[,k])
        nVec[eId]=nrow(tbl)
        if (exptIdVec[eId]==exptMultId) {
            nsnpThis=length(strsplit(tbl$genotype[1],split=" ")[[1]])
            nThis=nrow(tbl)
            genmatThis=matrix(NA,nThis,nsnpThis)
            for (i in 1:nThis) {genmatThis[i,]=strsplit(tbl$genotype[i],split=" ")[[1]]}
            genvecThis=rep(NA,nThis); names(genvecThis)=tbl$id
            for(i in 1:nThis) {genvecThis[i]=paste0(genmatThis[i,],collapse="")}
        }
    }

    if (truthCallType=="manual") {
        hugenList=NULL
    } else {
        ## Unique genotype combinations from singlet experiments
        hugenList=getUniqueGenoComb(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)
        hugenList=hugenList[which(names(hugenList)%in%paste0("SNACS",exptSingleId))]
    }
    
    genvec5=genvecThis
    hugen5=hugenList

    amb5=getAmbiguousGenotype(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)

    ####################################################################
    ## SNACS6

    exptMultId=6
    exptSingleId=c(2,3,4)
    exptIdVec=c(exptMultId,exptSingleId)

    nVec=rep(NA,length(exptIdVec))
    for (eId in 1:length(exptIdVec)) {
        tbl=read.table(paste0(dirData,"annCell_SNACS",exptIdVec[eId],exptNameSuffix,"_SNACS",exptMultId,"snp.txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        for (k in c("genotype","genoUnimputed")) tbl[,k]=as.character(tbl[,k])
        nVec[eId]=nrow(tbl)
        if (exptIdVec[eId]==exptMultId) {
            nsnpThis=length(strsplit(tbl$genotype[1],split=" ")[[1]])
            nThis=nrow(tbl)
            genmatThis=matrix(NA,nThis,nsnpThis)
            for (i in 1:nThis) {genmatThis[i,]=strsplit(tbl$genotype[i],split=" ")[[1]]}
            genvecThis=rep(NA,nThis); names(genvecThis)=tbl$id
            for(i in 1:nThis) {genvecThis[i]=paste0(genmatThis[i,],collapse="")}
        }
    }

    if (truthCallType=="manual") {
        hugenList=NULL
    } else {
        ## Unique genotype combinations from singlet experiments
        hugenList=getUniqueGenoComb(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)
        hugenList=hugenList[which(names(hugenList)%in%paste0("SNACS",exptSingleId))]
    }
    
    genvec6=genvecThis
    hugen6=hugenList

    amb6=getAmbiguousGenotype(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)

    ####################################################################
    ## SNACS7

    exptMultId=7
    exptSingleId=c(1,2,3,4)
    exptIdVec=c(exptMultId,exptSingleId)

    nVec=rep(NA,length(exptIdVec))
    for (eId in 1:length(exptIdVec)) {
        tbl=read.table(paste0(dirData,"annCell_SNACS",exptIdVec[eId],exptNameSuffix,"_SNACS",exptMultId,"snp.txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        for (k in c("genotype","genoUnimputed")) tbl[,k]=as.character(tbl[,k])
        nVec[eId]=nrow(tbl)
        if (exptIdVec[eId]==exptMultId) {
            nsnpThis=length(strsplit(tbl$genotype[1],split=" ")[[1]])
            nThis=nrow(tbl)
            genmatThis=matrix(NA,nThis,nsnpThis)
            for (i in 1:nThis) {genmatThis[i,]=strsplit(tbl$genotype[i],split=" ")[[1]]}
            genvecThis=rep(NA,nThis); names(genvecThis)=tbl$id
            for(i in 1:nThis) {genvecThis[i]=paste0(genmatThis[i,],collapse="")}
        }
    }

    if (truthCallType=="manual") {
        hugenList=NULL
    } else {
        ## Unique genotype combinations from singlet experiments
        hugenList=getUniqueGenoComb(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)
        hugenList=hugenList[which(names(hugenList)%in%paste0("SNACS",exptSingleId))]
    }
    
    genvec7=genvecThis
    hugen7=hugenList

    amb7=getAmbiguousGenotype(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)

    
}

    ####################################################################
    ## SNACS14

    exptMultId=14
    exptSingleId=c(105,106,127,128)
    exptIdVec=c(exptMultId,exptSingleId)

    nVec=rep(NA,length(exptIdVec))
    for (eId in 1:length(exptIdVec)) {
        tbl=read.table(paste0(dirData,"annCell_SNACS",exptIdVec[eId],exptNameSuffix,"_SNACS",exptMultId,"snp.txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
        for (k in c("genotype","genoUnimputed")) tbl[,k]=as.character(tbl[,k])
        nVec[eId]=nrow(tbl)
        if (exptIdVec[eId]==exptMultId) {
            nsnpThis=length(strsplit(tbl$genotype[1],split=" ")[[1]])
            nThis=nrow(tbl)
            genmatThis=matrix(NA,nThis,nsnpThis)
            for (i in 1:nThis) {genmatThis[i,]=strsplit(tbl$genotype[i],split=" ")[[1]]}
            genvecThis=rep(NA,nThis); names(genvecThis)=tbl$id
            for(i in 1:nThis) {genvecThis[i]=paste0(genmatThis[i,],collapse="")}
        }
    }

    if (truthCallType=="manual") {
        hugenList=NULL
    } else {
        ## Unique genotype combinations from singlet experiments
        hugenList=getUniqueGenoComb(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)
        hugenList=hugenList[which(names(hugenList)%in%paste0("SNACS",exptSingleId))]
    }
    
    genvec14=genvecThis
    hugen14=hugenList

    amb14=getAmbiguousGenotype(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)

    ####################################################################
    ## SNACS15

    if (maxSam>=8) {
        exptMultId=15
        exptSingleId=c(1,2,3,4,105,106,127,128)
        exptIdVec=c(exptMultId,exptSingleId)

        nVec=rep(NA,length(exptIdVec))
        for (eId in 1:length(exptIdVec)) {
            tbl=read.table(paste0(dirData,"annCell_SNACS",exptIdVec[eId],exptNameSuffix,"_SNACS",exptMultId,"snp.txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            for (k in c("genotype","genoUnimputed")) tbl[,k]=as.character(tbl[,k])
            nVec[eId]=nrow(tbl)
            if (exptIdVec[eId]==exptMultId) {
                nsnpThis=length(strsplit(tbl$genotype[1],split=" ")[[1]])
                nThis=nrow(tbl)
                genmatThis=matrix(NA,nThis,nsnpThis)
                for (i in 1:nThis) {genmatThis[i,]=strsplit(tbl$genotype[i],split=" ")[[1]]}
                genvecThis=rep(NA,nThis); names(genvecThis)=tbl$id
                for(i in 1:nThis) {genvecThis[i]=paste0(genmatThis[i,],collapse="")}
            }
        }

        if (truthCallType=="manual") {
            hugenList=NULL
        } else {
            ## Unique genotype combinations from singlet experiments
            hugenList=getUniqueGenoComb(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)
            hugenList=hugenList[which(names(hugenList)%in%paste0("SNACS",exptSingleId))]
        }
        
        genvec15=genvecThis
        hugen15=hugenList

        amb15=getAmbiguousGenotype(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)

    }

    ####################################################################
    ## SNACS16

    if (maxSam>=8) {
        exptMultId=16
        exptSingleId=c(1,2,3,4,105,106,127,128)
        exptIdVec=c(exptMultId,exptSingleId)

        nVec=rep(NA,length(exptIdVec))
        for (eId in 1:length(exptIdVec)) {
            tbl=read.table(paste0(dirData,"annCell_SNACS",exptIdVec[eId],exptNameSuffix,"_SNACS",exptMultId,"snp.txt"),sep="\t",h=T,quote="",comment.char="",as.is=T,fill=T)
            for (k in c("genotype","genoUnimputed")) tbl[,k]=as.character(tbl[,k])
            nVec[eId]=nrow(tbl)
            if (exptIdVec[eId]==exptMultId) {
                nsnpThis=length(strsplit(tbl$genotype[1],split=" ")[[1]])
                nThis=nrow(tbl)
                genmatThis=matrix(NA,nThis,nsnpThis)
                for (i in 1:nThis) {genmatThis[i,]=strsplit(tbl$genotype[i],split=" ")[[1]]}
                genvecThis=rep(NA,nThis); names(genvecThis)=tbl$id
                for(i in 1:nThis) {genvecThis[i]=paste0(genmatThis[i,],collapse="")}
            }
        }

        if (truthCallType=="manual") {
            hugenList=NULL
        } else {
            ## Unique genotype combinations from singlet experiments
            hugenList=getUniqueGenoComb(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)
            hugenList=hugenList[which(names(hugenList)%in%paste0("SNACS",exptSingleId))]
        }
        
        genvec16=genvecThis
        hugen16=hugenList

        amb16=getAmbiguousGenotype(exptName=paste0("SNACS",exptMultId),thresGenoCall=thresGenoCall,thresGenoNA=thresGenoNA)

    }

    ####################################################################
    ## Truth calls

    combos.2=function(single.list) {
        one=single.list[[1]]
        two=single.list[[2]]
        n1=length(one)
        n2=length(two)
        output=rep(NA,n1*n2)
        counter=0
        for(i in 1:n1)
        {
            seq1=as.numeric(strsplit(one[i],split="")[[1]])
            for(j in 1:n2)
            {
                counter=counter+1
                seq2=as.numeric(strsplit(two[j],split="")[[1]])
                seqsum=seq1+seq2
                seqsum[seqsum>1]=1
                output[counter]=paste0(seqsum,collapse="")
            }
        }
        output=sort(unique(output))
        output
    }

    combos.3=function(single.list) {
        one=single.list[[1]]
        two=single.list[[2]]
        three=single.list[[3]]
        output=rep(NA,length(one)*length(two)+length(one)*length(three)+length(two)*length(three)+length(one)*length(two)*length(three))

        counter=0
        for(i in 1:length(one))
        {
            seq1=as.numeric(unlist(strsplit(one[i],split="")))
            for(j in 1:length(two))
            {
                seq2=as.numeric(unlist(strsplit(two[j],split="")))
                s=seq1+seq2
                s[s>1]=1
                counter=counter+1
                output[counter]=paste0(s,collapse="")
                for(k in 1:length(three))
                {
                    seq3=as.numeric(unlist(strsplit(three[k],split="")))
                    s=seq1+seq3
                    s[s>1]=1
                    counter=counter+1
                    output[counter]=paste0(s,collapse="")
                    s=seq2+seq3
                    s[s>1]=1
                    counter=counter+1
                    output[counter]=paste0(s,collapse="")
                    s=seq1+seq2+seq3
                    s[s>1]=1
                    counter=counter+1
                    output[counter]=paste0(s,collapse="")
                }
            }
        }
        output=sort(unique(output))
        output
    }
    

    combos.4=function(single.list,fName=NULL) {
        #single.list=hugen14; fName=paste0(dirThis,"combGeno_",exptName,exptNameSuffix,"_",exptName,"snp")
        
        comb2Seqs=function(seqThis) {
            y=seqRest+seqThis
            y[y>1]=1
            apply(y,2,paste,collapse="")
        }

        single.list[[4]]=setdiff(single.list[[4]],single.list[[3]])
        single.list[[3]]=setdiff(single.list[[3]],single.list[[2]])


        n=unlist(lapply(single.list,length)); s=order(n)
        mat1=sapply(single.list[[s[1]]],function(x) {as.integer(strsplit(x,"")[[1]])},USE.NAMES=F)
        mat2=sapply(single.list[[s[2]]],function(x) {as.integer(strsplit(x,"")[[1]])},USE.NAMES=F)
        mat3=sapply(single.list[[s[3]]],function(x) {as.integer(strsplit(x,"")[[1]])},USE.NAMES=F)
        mat4=sapply(single.list[[s[4]]],function(x) {as.integer(strsplit(x,"")[[1]])},USE.NAMES=F)

        cat("\n\nNo. of geno: ",paste0(n,collapse=", "),"\n",sep="")

        output=NULL

        cat("\n\nOne N ",ncol(mat1),"\n",sep="")
        timeStamp3=Sys.time()
        for(i1 in 1:ncol(mat1)) {
        #for(i1 in 1) {
            cat(i1,"\n",sep="")

            output=NULL

            seq1=mat1[,i1]
            for(i2 in 1:ncol(mat2)) {
                seq2=mat2[,i2]
                s=seq1+seq2
                s[s>1]=1
                output=c(output,paste0(s,collapse=""))
                for(i3 in 1:ncol(mat3)) {
                    seq3=mat3[,i3]
                    
                    seqRest=matrix(c(seq1,seq2,seq1+seq2),ncol=3,byrow=F)
                    out=apply(mat3,2,comb2Seqs)
                    output=c(output,c(out))
                    
                    seqRest=matrix(c(seq1,seq2,seq3,seq1+seq2,seq1+seq3,seq2+seq3,seq1+seq2+seq3),ncol=7,byrow=F)
                    out=apply(mat4,2,comb2Seqs)
                    output=c(output,c(out))

                    out=apply(mat4,2,comb2Seqs)
                    
                    output=c(output,c(out))
                }
            }

            write.table(paste0("X",unique(output)),file=paste0(fName,"_",i1,".txt"),col.names=F,row.names=F, sep="\t",quote=F)

        }
        timeStamp3=c(timeStamp3,Sys.time())
        print(timeStamp3)
        print(diff(timeStamp3))

        #output=sort(unique(output))
        #output=sort(output)
        #output
    }

    combos.8=function(single.list) {
        #single.list=hugen15
        
        comb2Seqs=function(seqThis) {
            y=seqRest+seqThis
            y[y>1]=1
            apply(y,2,paste,collapse="")
        }
        
        #comb2Seqs=function(seqThis) {"a"}

        for (k in 3:length(single.list)) {single.list[[k]]=setdiff(single.list[[k]],single.list[[k-1]])}

        output=NULL

        matS=list()
        n=unlist(lapply(single.list,length)); s=order(n)
        for (k in 1:length(s)) {
            matS[[k]]=sapply(single.list[[s[k]]],function(x) {as.integer(strsplit(x,"")[[1]])},USE.NAMES=F)
        }

        cat("\n\nNo. of geno: ",paste0(n,collapse=", "),"\n",sep="")

        cat("\n\nOne N ",ncol(matS[[1]]),"\n")
        timeStamp3=Sys.time()
        for(i1 in 1:ncol(matS[[1]])) {
        #for(i1 in 1:5) {
            cat(i1,"\n")
            seq1=matS[[1]][,i1]
            for(i2 in 1:ncol(matS[[2]])) {
                seq2=matS[[2]][,i2]
                s=seq1+seq2
                s[s>1]=1
                output=c(output,paste0(s,collapse=""))
                for(i3 in 1:ncol(matS[[3]])) {
                    seq3=matS[[3]][,i3]
                    
                    seqRest=matrix(c(seq1,seq2,seq1+seq2),ncol=3,byrow=F)
                    out=apply(matS[[3]],2,comb2Seqs)
                    output=c(output,c(out))
                    
                    for(i4 in 1:ncol(matS[[4]])) {
                        seq4=matS[[4]][,i4]
                        
                        seqRest=matrix(c(seq1,seq2,seq3,seq1+seq2,seq1+seq3,seq2+seq3,seq1+seq2+seq3),ncol=7,byrow=F)
                        out=apply(matS[[4]],2,comb2Seqs)
                        output=c(output,c(out))
                        
                        for(i5 in 1:ncol(matS[[5]])) {
                            seq5=matS[[5]][,i5]
                            
                            seqRest=matrix(c(seq1,seq2,seq3,seq4,
                                seq1+seq2,seq1+seq3,seq1+seq4,seq2+seq3,seq2+seq4,seq3+seq4,
                                seq1+seq2+seq3,seq1+seq2+seq4,seq1+seq2+seq3+seq4),ncol=13,byrow=F)
                            out=apply(matS[[5]],2,comb2Seqs)
                            output=c(output,c(out))
                            
                            for(i6 in 1:ncol(matS[[6]])) {
                                seq6=matS[[6]][,i6]
                                
                                seqRest=matrix(c(seq1,seq2,seq3,seq4,seq5,
                                    seq1+seq2,seq1+seq3,seq1+seq4,seq1+seq5,
                                    seq2+seq3,seq2+seq4,seq2+seq5,
                                    seq3+seq4,seq3+seq5,seq4+seq5,
                                    seq1+seq2+seq3,seq1+seq2+seq4,seq1+seq2+seq5,
                                    seq2+seq3+seq4,seq2+seq3+seq5,seq3+seq4+seq5,
                                    seq1+seq2+seq3+seq4,seq1+seq2+seq3+seq5,
                                    seq2+seq3+seq4+seq5,seq1+seq2+seq3+seq4+seq5),ncol=25,byrow=F)
                                out=apply(matS[[6]],2,comb2Seqs)
                                output=c(output,c(out))
                                
                                for(i7 in 1:ncol(matS[[7]])) {
                                    seq7=matS[[7]][,i7]
                                    
                                    seqRest=matrix(c(seq1,seq2,seq3,seq4,seq5,seq6,
                                        seq1+seq2,seq1+seq3,seq1+seq4,seq1+seq5,seq1+seq6,
                                        seq2+seq3,seq2+seq4,seq2+seq5,seq2+seq6,
                                        seq3+seq4,seq3+seq5,seq3+seq6,seq4+seq5,seq4+seq6,seq5+seq6,
                                        seq1+seq2+seq3,seq1+seq2+seq4,seq1+seq2+seq5,seq1+seq2+seq6,
                                        seq2+seq3+seq4,seq2+seq3+seq5,seq2+seq3+seq6,
                                        seq3+seq4+seq5,seq3+seq4+seq6,seq4+seq5+seq6,
                                        seq1+seq2+seq3+seq4,seq1+seq2+seq3+seq5,seq1+seq2+seq3+seq6,
                                        seq2+seq3+seq4+seq5,seq2+seq3+seq4+seq6,seq3+seq4+seq5+seq6,
                                        seq1+seq2+seq3+seq4+seq5,seq1+seq2+seq3+seq4+seq6,seq2+seq3+seq4+seq5+seq6,
                                        seq1+seq2+seq3+seq4+seq5+seq6),ncol=41,byrow=F)
                                    out=apply(matS[[7]],2,comb2Seqs)
                                    output=c(output,c(out))
                                    
                                    seqRest=matrix(c(seq1,seq2,seq3,seq4,seq5,seq6,seq7,
                                        seq1+seq2,seq1+seq3,seq1+seq4,seq1+seq5,seq1+seq6,seq1+seq7,
                                        seq2+seq3,seq2+seq4,seq2+seq5,seq2+seq6,seq2+seq7,
                                        seq3+seq4,seq3+seq5,seq3+seq6,seq3+seq7,
                                        seq4+seq5,seq4+seq6,seq4+seq7,
                                        seq5+seq6,seq5+seq7,seq6+seq7,
                                        seq1+seq2+seq3,seq1+seq2+seq4,seq1+seq2+seq5,seq1+seq2+seq6,seq1+seq2+seq7,
                                        seq2+seq3+seq4,seq2+seq3+seq5,seq2+seq3+seq6,seq2+seq3+seq7,
                                        seq3+seq4+seq5,seq3+seq4+seq6,seq3+seq4+seq7,
                                        seq4+seq5+seq6,seq4+seq5+seq7,seq5+seq6+seq7,
                                        seq1+seq2+seq3+seq4,seq1+seq2+seq3+seq5,seq1+seq2+seq3+seq6,seq1+seq2+seq3+seq7,
                                        seq2+seq3+seq4+seq5,seq2+seq3+seq4+seq6,seq2+seq3+seq4+seq7,
                                        seq3+seq4+seq5+seq6,seq3+seq4+seq5+seq7,seq4+seq5+seq6+seq7,
                                        seq1+seq2+seq3+seq4+seq5,seq1+seq2+seq3+seq4+seq6,seq1+seq2+seq3+seq4+seq7,
                                        seq2+seq3+seq4+seq5+seq6,seq2+seq3+seq4+seq5+seq7,
                                        seq3+seq4+seq5+seq6+seq7,seq1+seq2+seq3+seq4+seq5+seq6,seq1+seq2+seq3+seq4+seq5+seq7,
                                        seq2+seq3+seq4+seq5+seq6+seq7,seq1+seq2+seq3+seq4+seq5+seq6+seq7),ncol=63,byrow=F)
                                    out=apply(matS[[8]],2,comb2Seqs)
                                    output=c(output,c(out))
                                }
                            }
                        }
                    }
                }
            }
        }
        
        timeStamp3=c(timeStamp3,Sys.time())
        print(timeStamp3)
        print(diff(timeStamp3))

        output=sort(unique(output))
        output

    }

    ###Truth function based on genotypes matching single experiment data to be used estimate accuracy
    ###Requires no missing genotypic data, either due to imputation or filtering
    ###gv is the vector of genotypes of the multiplexed samples as a single string
    ###single.list is a list of genotypes from the single experiment as a single string
    ###output is predicted sample
    ###genotypes look like string with no spaces: 011100100
    ###not: 0 1 0 1 1 1 0 0
    ###to turn string with spaces into one with no spaces run paste0(string,collapse="")
    ###match.combos.extended are to find multiplets defined as genotype with mutations no seen in an individual sample
    truth.function=function(gv,single.list,amb=NULL,exptName="",loadFlag=F) {
        #gv=genvec6; single.list=list(hugenvec26,hugenvec36,hugenvec46); amb=amb6
        #gv=genvec7; single.list=list(hugenvec17,hugenvec27,hugenvec37,hugenvec47); amb=amb7
        #gv=genvec14; single.list=hugen14; amb=amb14; exptName="SNACS14"
        
        cat("\n------------ truth.function: ",exptName," --------------\n",sep="")
        timeStamp=Sys.time()
        
    ###How many samples have been multiplexed
        n.samples=length(single.list)
        if (loadFlag) {
            dirThis="../output/accuracy/combGeno/"
            fileList=dir(dirThis,pattern=paste0("combGeno_",exptName,exptNameSuffix,"_",exptName,"snp"))
            combos=NULL
            for (fId in 1:length(fileList)) {
                dat=read.table(paste0(dirThis,fileList[fId]),sep="\t",h=F,quote="",comment.char="",as.is=T,fill=T)
                combos=c(combos,dat[,1])
            }
            combos=sub("X","",sort(unique(combos)))
        } else {
            if(n.samples==2) combos=combos.2(single.list)
            if(n.samples==3) combos=combos.3(single.list)
            if(n.samples==4) combos=combos.4(single.list)
            if(n.samples==8) combos=combos.8(single.list)
        }
        n.combos=length(combos)
    ###ugv are the unique genotypes
        ugv=sort(unique(gv))
        n.genotypes=length(ugv)
    ###Unlist the singles for matching purposes
        singles=unlist(single.list)
        n.singles=length(singles)
        match.singles.exact=match(ugv,singles)
        match.combos.exact=match(ugv,combos)
        no.matches=which(is.na(match.singles.exact) & is.na(match.combos.exact))
    ###Transform the strings to be separate components to make comparisons
        singles.matrix=matrix(NA,n.singles,length(unlist(strsplit(singles[1],split=""))))
        for(i in 1:length(singles)) singles.matrix[i,]=unlist(strsplit(singles[i],split=""))
        #combos.matrix=matrix(NA,length(combos),length(unlist(strsplit(combos[1],split=""))))
        #for(i in 1:n.combos) combos.matrix[i,]=unlist(strsplit(combos[i],split=""))
    ###Run through the genotypes to see if they can be multplets
    ###They are multiplets if somewhere they have a 1 where a single has a 0 for every single
        match.combos.extended=rep(FALSE,n.genotypes)
        for(i in 1:length(no.matches))
        {
            match.combos.extended.current=rep(FALSE,n.singles)
            new.i=no.matches[i]
            new.genotype=unlist(strsplit(ugv[new.i],split=""))
            for(j in 1:n.singles)
            {
                new.single=singles.matrix[j,]
                if(length(which(new.genotype==1 & new.single==0)>0)) match.combos.extended.current[j]=TRUE
            }
            if(sum(match.combos.extended.current)==n.singles) match.combos.extended[new.i]=TRUE
        }
    ###Combine exact matches with extended matches for multiplets
        indicator.multiplet=rep(FALSE,n.genotypes)
        indicator.multiplet[which(!is.na(match.combos.exact)|match.combos.extended)]=TRUE
    ## A true ambiguous call for a cell is when the cell is not a multiplet but cannot be uniquely assigned to a single sample
        match.ambiguous.exact=rep(FALSE,n.genotypes)
        match.ambiguous.exact[which(ugv%in%singles[duplicated(singles)])]=TRUE
        #match.ambiguous.exact[which(indicator.multiplet)]=FALSE

        ambThis=amb
        
        ambThis$genotype=gsub(" ","",amb$genotype)
    ###Now look for singles among those not multiplet or already single, commented out for now
    ###    no.single.multiplet=which(indicator.multiplet==FALSE & is.na(match.singles.exact))
    ###    print(no.single.multiplet)
    ###    match.singles.extended=rep(FALSE,n.genotypes)
    ###    if(length(no.single.multiplet)>0)
    ###    {
    ###        for(i in 1:length(no.single.multiplet))
    ###        {
    ###            match.singles.extended.current=rep(FALSE,n.singles)
    ###            new.i=no.single.multiplet[i]
    ###            new.genotype=unlist(strsplit(ugv[new.i],split=""))
    ###            for(j in 1:n.singles)
    ###            {
    ###                new.single=singles.matrix[j,]
    ###                if(length(which(new.genotype==1 & new.single==0)==0)) match.singles.extended.current[j]=TRUE
    ###            }
    ###            if(sum(match.singles.extended.current)>0) match.singles.extended[new.i]=TRUE
    ###        }
    ###    }
    ###Combine exact matches with extended matches for singles
    ###    indicator.single=rep(FALSE,n.genotypes)
    ###    indicator.single[which(!is.na(match.singles.exact)|match.singles.extended)]=TRUE
    ###Now make final calls based on unique genotypes
        singles.patterns=rep(NA,n.samples)
        for(i in 1:n.samples) singles.patterns[i]=length(single.list[[i]])
        singles.ends=cumsum(singles.patterns)
        singles.starts=c(1,singles.ends[-n.samples]+1)
        match.gv.ugv=match(gv,ugv)
        matches.vector=rep(NA,n.genotypes)
        matches.vector[indicator.multiplet]="multiplet"
        for(i in 1:n.samples)
        {
            ii=which(match.singles.exact>=singles.starts[i] & match.singles.exact<=singles.ends[i])
            #ii=which(match.singles.exact>=singles.starts[i] & match.singles.exact<=singles.ends[i] & !indicator.multiplet)
            nm=paste0("Sample.",i)
            matches.vector[ii]=nm
            if (!is.null(amb)) {
                j=singles.starts[i]:singles.ends[i]
                #jj=which(amb$genotype%in%singles[j] & amb$genotype%in%singles[-j] & amb$samId1==nm)
                jj=which(amb$genotype%in%singles[j] & (amb$samId1==nm | amb$samId2==nm))
                if (length(jj)!=0) {
                    i1=which(ugv[ii]%in%amb$genotype[jj])
                    if (length(i1)!=0) matches.vector[ii[i1]]="ambiguous"
                }
            }
        }
        matches.vector[match.ambiguous.exact]="ambiguous"
        final.calls=matches.vector[match.gv.ugv]
        #output=list(ugv=ugv,match.singles.exact=match.singles.exact,match.combos.exact=match.combos.exact,no.matches=no.matches,singles=singles,singles.matrix=singles.matrix,combos.matrix=combos.matrix,match.combos.extended=match.combos.extended,indicator.multiplet=indicator.multiplet,match.gv.ugv=match.gv.ugv,matches.vector=matches.vector,final.calls=final.calls)

        timeStamp=c(timeStamp,Sys.time())
        print(timeStamp)
        print(diff(timeStamp))

        #output
        final.calls
    }
    
    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))

    dirThis="../output/accuracy/combGeno/"
    if (!file.exists(dirThis)) dir.create(file.path(dirThis))

    if (F) {
    ## ------------------------------
    exptName="SNACS5"
    
    cat("\n\n------------ ",exptName,sep="")

    if (saveFlag) {
        tbl=paste0("X",combos.2(hugen5))
        write.table(tbl,file=paste0(dirThis,"combGeno_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=F,row.names=F, sep="\t",quote=F)
    }

    foo5=truth.function(genvec5,hugen5,amb=amb5,exptName=exptName,loadFlag=loadFlag)
    
    #makeHashPairScatterPlot(exptName=exptName,exptNameSuffix)
    tbl=data.frame(id=names(genvec5))
    tbl$hashCall_truth=foo5
    write.table(tbl,file=paste0(dirOutput,"trueHashCall_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=T,row.names=F, sep="\t",quote=F)

    ## ------------------------------
    exptName="SNACS6"
    
    cat("\n\n------------ ",exptName,sep="")

    if (saveFlag) {
        tbl=paste0("X",combos.3(hugen6))
        write.table(tbl,file=paste0(dirThis,"combGeno_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=F,row.names=F, sep="\t",quote=F)
    }

    foo6=truth.function(genvec6,hugen6,amb=amb6,exptName=exptName,loadFlag=loadFlag)
    
    tbl=data.frame(id=names(genvec6))
    tbl$hashCall_truth=foo6
    write.table(tbl,file=paste0(dirOutput,"trueHashCall_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=T,row.names=F, sep="\t",quote=F)

    ## ------------------------------
    exptName="SNACS7"
    
    cat("\n\n------------ ",exptName,sep="")

    if (saveFlag) {
        tbl=paste0("X",combos.4(hugen7))
        write.table(tbl,file=paste0(dirThis,"combGeno_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=F,row.names=F, sep="\t",quote=F)
    }

    foo7=truth.function(genvec7,hugen7,amb=amb7,exptName=exptName,loadFlag=loadFlag)

    tbl=data.frame(id=names(genvec7))
    tbl$hashCall_truth=foo7
    write.table(tbl,file=paste0(dirOutput,"trueHashCall_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=T,row.names=F, sep="\t",quote=F)

    }

    ## ------------------------------
    if (T) {

    exptName="SNACS14"
    
    cat("\n\n------------ ",exptName,sep="")

    if (saveFlag) {
        combos.4(hugen14,fName=paste0(dirThis,"combGeno_",exptName,exptNameSuffix,"_",exptName,"snp"))
        #tbl=paste0("X",combos.4(hugen14))
        #write.table(tbl,file=paste0(dirThis,"combGeno_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=F,row.names=F, sep="\t",quote=F)
    }


    foo14=truth.function(genvec14,hugen14,amb=amb14,exptName=exptName,loadFlag=loadFlag)

    tbl=data.frame(id=names(genvec14))
    tbl$hashCall_truth=foo14
    write.table(tbl,file=paste0(dirOutput,"trueHashCall_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=T,row.names=F, sep="\t",quote=F)

    }
    ## ------------------------------
    if (maxSam>=8) {
        if (F) {


    exptName="SNACS15"
    
    cat("\n\n------------ ",exptName,sep="")

    if (saveFlag) {
        tbl=paste0("X",combos.8(hugen15))
        write.table(tbl,file=paste0(dirThis,"combGeno_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=F,row.names=F, sep="\t",quote=F)
    }

    foo15=truth.function(genvec15,hugen15,amb=amb15,exptName=exptName,loadFlag=loadFlag)

    tbl=data.frame(id=names(genvec15))
    tbl$hashCall_truth=foo15
    write.table(tbl,file=paste0(dirOutput,"trueHashCall_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=T,row.names=F, sep="\t",quote=F)
    
    }
    ## ------------------------------
    
    exptName="SNACS16"

    cat("\n\n------------ ",exptName,sep="")

    if (saveFlag) {
        tbl=paste0("X",combos.8(hugen16))
        write.table(tbl,file=paste0(dirThis,"combGeno_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=F,row.names=F, sep="\t",quote=F)
    }

    foo16=truth.function(genvec16,hugen16,amb=amb16,exptName=exptName,loadFlag=loadFlag)

    tbl=data.frame(id=names(genvec16))
    tbl$hashCall_truth=foo16
    write.table(tbl,file=paste0(dirOutput,"trueHashCall_",exptName,exptNameSuffix,"_",exptName,"snp.txt"),col.names=T,row.names=F, sep="\t",quote=F)


    }
    ## ------------------------------

    timeStamp=Sys.time()
    print(timeStamp)

}

####################################################################
####################################################################
