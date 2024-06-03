## SNACS_Manuscript

####################################################################
####################################################################
## Get unique logical genotype-combinations from all genotype-combinations from an experiment
getUniqueGenoComb=function(exptName,thresGenoCall=0,dirData="../output/accuracy/snpPosPropMat/") {
    cat("\n------------ getUniqueGenoComb: ",exptName," --------------\n",sep="")
    timeStamp=Sys.time()
    
    fName=paste0(dirData,"snpPosPropMat_",exptName,exptNameSuffix,"_",exptName,"snp.txt")
    if (file.exists(fName)) {
        tbl=read.delim(fName,sep="\t",header=TRUE,as.is=TRUE)
        snpMat=as.matrix(tbl[nrow(tbl):1,4:ncol(tbl)])
        snpMat[snpMat<=thresGenoCall]=0
        snpMat[snpMat>=(1-thresGenoCall)]=1
        snpMat[snpMat>thresGenoCall & snpMat<(1-thresGenoCall)]=0.5
        genoUniq=list()
        for (k in 1:ncol(snpMat)) {
            x=snpMat[,k]
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

getAmbiguousGenotype=function(exptName,thresGenoCall=0,dirData="../output/accuracy/snpPosPropMat/") {
    cat("\n------------ getAmbiguousGenotype: ",exptName," --------------\n",sep="")
    timeStamp=Sys.time()
    
    fName=paste0(dirData,"snpPosPropMat_",exptName,exptNameSuffix,"_",exptName,"snp.txt")
    if (file.exists(fName)) {
        tbl=read.delim(fName,sep="\t",header=TRUE,as.is=TRUE)
        snpMat=as.matrix(tbl[nrow(tbl):1,4:ncol(tbl)])
        snpMat[snpMat<=thresGenoCall]=0
        snpMat[snpMat>=(1-thresGenoCall)]=1
        snpMat[snpMat>thresGenoCall & snpMat<(1-thresGenoCall)]=0.5
        genoUniq=list()
        for (k in 1:ncol(snpMat)) {
            x=snpMat[,k]
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
    if (exptName=="SNACS5") {x=cbind(dat5,foo);
    } else if (exptName=="SNACS6") {x=cbind(dat6,foo=foo2)
    } else if (exptName=="SNACS7") {x=cbind(dat7,foo=foo3)}
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

makeTruthCall=function(exptNameSuffix,thresGenoCall,dirData="../output/accuracy/annCell/",dirOutput="../output/accuracy/trueHashCall/") {
    #exptNameSuffix=exptNameSuffix; thresGenoCall=thresGenoCall; dirData="../output/accuracy/annCell/"; dirOutput="../output/accuracy/trueHashCall/"
    
    #truthCallType="manual" ## NOT used. For hand curated singlet truth calls
    truthCallType="automatic"
    
    timeStamp=Sys.time()
    print(timeStamp)

    cat("\n\n---------------------------------\n")
    cat("------------------",exptNameSuffix," thres ",thresGenoCall," ---------------\n",sep="")
    cat("---------------------------------\n\n")

    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))

    
    #dirOutput=paste0("../output/accuracy/truthCall",ifelse(truthCallType=="manual","Manual",paste0("Auto_thres",thresGenoCall)),"/trueHashCall/")

    ###
    ###
    ###SNACS5, done like SNACS6 and SNACS7

    dat5=read.delim(paste0(dirData,"annCell_SNACS5",exptNameSuffix,"_SNACS5snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    for (k in c("genotype","genoUnimputed")) dat5[,k]=as.character(dat5[,k])
    nsnp5=length(strsplit(dat5$genotype[1],split=" ")[[1]])
    n5=nrow(dat5)
    
    #amb5=getAmbiguousGenotype(exptName="SNACS5",hashNames=unique(dat5$snacsRnd2))
    amb5=getAmbiguousGenotype(exptName="SNACS5",thresGenoCall=thresGenoCall)

    ###SNACS1, SNACS2

    dat15=read.delim(paste0(dirData,"annCell_SNACS1",exptNameSuffix,"_SNACS5snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    dat25=read.delim(paste0(dirData,"annCell_SNACS2",exptNameSuffix,"_SNACS5snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    for (k in c("genotype","genoUnimputed")) dat15[,k]=as.character(dat15[,k])
    for (k in c("genotype","genoUnimputed")) dat25[,k]=as.character(dat25[,k])
    n15=nrow(dat15)
    n25=nrow(dat25)

    ###genmat6 is the genotype of SNACS6

    genmat5=matrix(NA,n5,nsnp5)

    for(i in 1:n5)
    {
        genmat5[i,]=strsplit(dat5$genotype[i],split=" ")[[1]]
    }

    ###genmat15, genmat35 and genmat45 are genotypes of components of genmat5

    genmat15=matrix(NA,n15,nsnp5)
    for(i in 1:n15)
    {
        genmat15[i,]=strsplit(dat15$genotype[i],split=" ")[[1]]
    }

    genmat25=matrix(NA,n25,nsnp5)
    for(i in 1:n25)
    {
        genmat25[i,]=strsplit(dat25$genotype[i],split=" ")[[1]]
    }

    ###genvec collapses the genotypes

    genvec5=rep(NA,n5)
    for(i in 1:n5)
    {
        genvec5[i]=paste0(genmat5[i,],collapse="")
    }

    genvec15=rep(NA,n15)
    for(i in 1:n15)
    {
        genvec15[i]=paste0(genmat15[i,],collapse="")
    }

    genvec25=rep(NA,n25)
    for(i in 1:n25)
    {
        genvec25[i]=paste0(genmat25[i,],collapse="")
    }

    ###Make ugenvec as unique components of genvec

    ugenvec5=sort(unique(genvec5))
    ugenvec15=sort(unique(genvec15))
    ugenvec25=sort(unique(genvec25))

    ###Do the components ever match?

    match.15.25=match(ugenvec15,ugenvec25)

    ###No, components are unique

    ###Count up the calls in SNACS5

    which.51=which(dat5$snacsRnd2=="TS.1")
    which.52=which(dat5$snacsRnd2=="TS.2")
    which.512=which(dat5$snacsRnd2=="TS.1_TS.2")

    ###1587 called 1, 765 called 2, 275 called 12

    genvec.51=genvec5[which.51]
    genvec.52=genvec5[which.52]
    genvec.512=genvec5[which.512]
    genvec.5multiplet=genvec.512

    ###How many match expected genotypes for patient 1

    match.51.15=match(genvec.51,ugenvec15)

    ###3 of 1587 (0.2%) called patient 1 don't match

    match.51.25=match(genvec.51,ugenvec25)

    ###0 out of 1587 (0%) matches patient 2

    ###How many match expected genotypes for patient 2

    match.52.25=match(genvec.52,ugenvec25)

    ###23 of 765 (3.0%) called patient 2 don't match

    match.52.15=match(genvec.52,ugenvec15)

    ###0 out of 765 (0%) matches patient 1

    ###Now which multiplets match 1, 2

    match.5multiplet.15=match(genvec.5multiplet,ugenvec15)
    match.5multiplet.25=match(genvec.5multiplet,ugenvec25)

    ###6 of 275 (2.2%) matches patient 2, 21 of 275 (7.6%)

    ###Logical ordering of genotypes

    ###Patient 1, 0EEEEE10, so 00000010, 00000110,
    ###Patient 2, 1000000E, so 10000000, 10000001


    if (truthCallType=="manual") {
        ## Not used
        ## Related to truth call based on 2 rounds of SNP clustering
        ###Adding an "h" for hand curated genotypes
        hugenvec15=c("10000000","10000001")
        hugenvec25=c("00000010","00000110","00001010","00001110","00010010","00010110","00011010","00011110","00100010","01000110","01010010","00101110","00110010","00110110","00111010","00111110","01000010","01000110","01001010","01001110","01010010","01010110","01011010","01011110","01100010","01100110","01101010","01101110","01110010","01110110","01111010","01111110")
    } else {
        ## Unique genotype combinations from singlet experiments
        res=getUniqueGenoComb(exptName="SNACS5",thresGenoCall=thresGenoCall)
        hugenvec15=res[["SNACS1"]]
        hugenvec25=res[["SNACS2"]]
    }

    ###Now match to patient 1

    hmatch.51.15=match(genvec.51,hugenvec15)
    hmatch.51.25=match(genvec.51,hugenvec25)

    ###55 of 885 (5.3%) don't match, 0 (0%) matches patient 3 and 1 (0.1%) matches patient 4

    ###Now match to patient 2

    hmatch.52.15=match(genvec.52,hugenvec15)
    hmatch.52.25=match(genvec.52,hugenvec25)

    ###3 of 574 (0.5%) don't match, 0 (0%) matches patient 2 and 0 (0%) matches patient 4

    ###Now which multiplets match 2, 3 and 4

    hmatch.5multiplet.15=match(genvec.5multiplet,hugenvec15)
    hmatch.5multiplet.25=match(genvec.5multiplet,hugenvec25)

    ###13 (1.5%) match patient 2, 10 (1.3%) match patient 3, and 25 (4.4%) match patient 4
    ###Check whether consistent with multiplets of other patients

    ###Find codes of multiplets that matched individual patients, later check if they match combos

    hugenvec.5multiplet.1=hugenvec15[hmatch.5multiplet.15[!is.na(hmatch.5multiplet.15)]]
    hugenvec.5multiplet.2=hugenvec25[hmatch.5multiplet.25[!is.na(hmatch.5multiplet.25)]]

    ###Now find all the combinations

    h5.combos=rep(NA,length(hugenvec15)*length(hugenvec25))

    count=0
    for(i in 1:length(hugenvec15))
    {
        one=as.numeric(unlist(strsplit(hugenvec15[i],split="")))
        for(j in 1:length(hugenvec25))
            {
                two=as.numeric(unlist(strsplit(hugenvec25[j],split="")))
                s=one+two
                s[s>1]=1
                count=count+1
                h5.combos[count]=paste0(s,collapse="")
            }
    }

    h5.combos=sort(unique(h5.combos))

    ###match multiplets to combos

    hmatch.5multiplet.combos=match(genvec.5multiplet,h5.combos)

    ###59 out of 275 don't match (21.5%).  There are 68 genotypes from multiplets

    genvec.5multiplet.nomatch=genvec.5multiplet[is.na(hmatch.5multiplet.combos)]

    ###Do multiplets that match individual patients also match combos

    match.hugenvec.5multiplet.1.combos=match(hugenvec.5multiplet.1,h5.combos)
    match.hugenvec.5multiplet.2.combos=match(hugenvec.5multiplet.2,h5.combos)

    ###13 of 13 multiplets that match patient 2 also match combos, 0 of 10 multiplets that match patient 3 match combos, and 0 of 25 that match patient 4 match combos

    ###
    ###
    ###SNACS6

    dat6=read.delim(paste0(dirData,"annCell_SNACS6",exptNameSuffix,"_SNACS6snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    for (k in c("genotype","genoUnimputed")) dat6[,k]=as.character(dat6[,k])
    nsnp6=length(strsplit(dat6$genotype[1],split=" ")[[1]])
    n6=nrow(dat6)
    
    #amb6=getAmbiguousGenotype(exptName="SNACS6",hashNames=unique(dat6$snacsRnd2))
    amb6=getAmbiguousGenotype(exptName="SNACS6",thresGenoCall=thresGenoCall)


    ###SNACS2, SNACS3, SNACS4

    dat26=read.delim(paste0(dirData,"annCell_SNACS2",exptNameSuffix,"_SNACS6snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    dat36=read.delim(paste0(dirData,"annCell_SNACS3",exptNameSuffix,"_SNACS6snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    dat46=read.delim(paste0(dirData,"annCell_SNACS4",exptNameSuffix,"_SNACS6snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    for (k in c("genotype","genoUnimputed")) dat26[,k]=as.character(dat26[,k])
    for (k in c("genotype","genoUnimputed")) dat36[,k]=as.character(dat36[,k])
    for (k in c("genotype","genoUnimputed")) dat46[,k]=as.character(dat46[,k])
    n26=nrow(dat26)
    n36=nrow(dat36)
    n46=nrow(dat46)

    ###genmat6 is the genotype of SNACS6

    genmat6=matrix(NA,n6,nsnp6)

    for(i in 1:n6)
    {
        genmat6[i,]=strsplit(dat6$genotype[i],split=" ")[[1]]
    }

    ###genmat26, genmat36 and genmat46 are genotypes of components of genmat6

    genmat26=matrix(NA,n26,nsnp6)
    for(i in 1:n26)
    {
        genmat26[i,]=strsplit(dat26$genotype[i],split=" ")[[1]]
    }

    genmat36=matrix(NA,n36,nsnp6)
    for(i in 1:n36)
    {
        genmat36[i,]=strsplit(dat36$genotype[i],split=" ")[[1]]
    }

    genmat46=matrix(NA,n46,nsnp6)
    for(i in 1:n46)
    {
        genmat46[i,]=strsplit(dat46$genotype[i],split=" ")[[1]]
    }

    ###genvec collapses the genotypes

    genvec6=rep(NA,n6)
    for(i in 1:n6)
    {
        genvec6[i]=paste0(genmat6[i,],collapse="")
    }

    genvec26=rep(NA,n26)
    for(i in 1:n26)
    {
        genvec26[i]=paste0(genmat26[i,],collapse="")
    }

    genvec36=rep(NA,n36)
    for(i in 1:n36)
    {
        genvec36[i]=paste0(genmat36[i,],collapse="")
    }

    genvec46=rep(NA,n46)
    for(i in 1:n46)
    {
        genvec46[i]=paste0(genmat46[i,],collapse="")
    }

    ###Make ugenvec as unique components of genvec

    ugenvec6=sort(unique(genvec6))
    ugenvec26=sort(unique(genvec26))
    ugenvec36=sort(unique(genvec36))
    ugenvec46=sort(unique(genvec46))

    ###Do the components ever match?

    match.26.36=match(ugenvec26,ugenvec36)
    match.26.46=match(ugenvec26,ugenvec46)
    match.36.46=match(ugenvec36,ugenvec46)

    ###No, components are unique

    ###Count up the calls in SNACS6

    which.62=which(dat6$snacsRnd2=="TS.2")
    which.63=which(dat6$snacsRnd2=="TS.3")
    which.64=which(dat6$snacsRnd2=="TS.4")
    which.623=which(dat6$snacsRnd2=="TS.2_TS.3")
    which.624=which(dat6$snacsRnd2=="TS.2_TS.4")
    which.634=which(dat6$snacsRnd2=="TS.3_TS.4")
    which.6234=which(dat6$snacsRnd2=="TS.2_TS.3_TS.4")

    ###885 called 2, 574 called 3, 1093 called 4, 196 called 23, 79 called 24, 191 called 34, 326 called 234

    genvec.62=genvec6[which.62]
    genvec.63=genvec6[which.63]
    genvec.64=genvec6[which.64]
    genvec.623=genvec6[which.623]
    genvec.624=genvec6[which.624]
    genvec.634=genvec6[which.634]
    genvec.6234=genvec6[which.6234]
    genvec.6multiplet=genvec6[c(which.623,which.624,which.634,which.6234)]

    ###How many match expected genotypes for patient 2

    match.62.26=match(genvec.62,ugenvec26)

    ###38 of 885 (4.3%) called patient 2 don't match

    match.62.36=match(genvec.62,ugenvec36)
    match.62.46=match(genvec.62,ugenvec46)

    ###0 out of 885 (0%) matches patient 3 and 1 out of 885 (.1%) matches patient 4

    ###How many match expected genotypes for patient 3

    match.63.36=match(genvec.63,ugenvec36)

    ###3 of 574 (0.5%) called patient 3 don't match

    match.63.26=match(genvec.63,ugenvec26)
    match.63.46=match(genvec.63,ugenvec46)

    ###0 out of 574 (0%) matches patient 2 and 0 out of 574 (0.1%) matches patient 4

    ###How many match expected genotypes for patient 4

    match.64.46=match(genvec.64,ugenvec46)

    ###13 of 1093 (1.2%) called patient 4 don't match

    match.64.26=match(genvec.64,ugenvec26)
    match.64.36=match(genvec.64,ugenvec36)

    ###3 out of 1093 (0.3%) matches patient 2 and 0 out of 1093 (0%) matches patient 3

    ###Now which multiplets match 2, 3 and 4

    match.6multiplet.26=match(genvec.6multiplet,ugenvec26)
    match.6multiplet.36=match(genvec.6multiplet,ugenvec36)
    match.6multiplet.46=match(genvec.6multiplet,ugenvec46)

    ###189 of 792 (23.4%) matches patient 2, 103 of 792 match patient 3 (13.0%), 73 of 792 match patient 3 (9.2%)

    ###Logical ordering of genotypes

    ###Patient 2, E11E001, so 0110001 (15), 0111001 (122), 1110001 (95), 1111001 (1461), all the rest (12)
    ###Patient 3, 0110010 (1841), all the rest (16)
    ###Patient 4, 0000E01, so 0000001 (143), 0000101 (5373), all the rest (33)

    if (truthCallType=="manual") {
        ## NOT USED
        ## Related to truth call based on 2 rounds of SNP clustering
        ###Adding an "h" for hand curated genotypes
        hugenvec26=c("0110001","0111001","1110001","1111001")
        hugenvec36="0110010"
        hugenvec46=c("0000001","0000101")
    } else {
        ## Unique genotype combinations from singlet experiments
        res=getUniqueGenoComb(exptName="SNACS6",thresGenoCall=thresGenoCall)
        hugenvec26=res[["SNACS2"]]
        hugenvec36=res[["SNACS3"]]
        hugenvec46=res[["SNACS4"]]
    }

    ###Now match to patient 2

    hmatch.62.26=match(genvec.62,hugenvec26)
    hmatch.62.36=match(genvec.62,hugenvec36)
    hmatch.62.46=match(genvec.62,hugenvec46)

    ###56 of 885 (6.3%) don't match, 0 (0%) matches patient 3 and 1 (0.1%) matches patient 4

    ###Now match to patient 3

    hmatch.63.26=match(genvec.63,hugenvec26)
    hmatch.63.36=match(genvec.63,hugenvec36)
    hmatch.63.46=match(genvec.63,hugenvec46)

    ###3 of 574 (0.5%) don't match, 0 (0%) matches patient 2 and 0 (0%) matches patient 4

    ###Now match to patient 4

    hmatch.64.26=match(genvec.64,hugenvec26)
    hmatch.64.36=match(genvec.64,hugenvec36)
    hmatch.64.46=match(genvec.64,hugenvec46)

    ###55 of 1093 (5.0%) don't match, 0 (0%) matches patient 2 and 0 (0%) matches patient 4

    ###Now which multiplets match 2, 3 and 4

    hmatch.6multiplet.26=match(genvec.6multiplet,hugenvec26)
    hmatch.6multiplet.36=match(genvec.6multiplet,hugenvec36)
    hmatch.6multiplet.46=match(genvec.6multiplet,hugenvec46)

    ###13 (1.6%) match patient 2, 10 (1.3%) match patient 3, and 35 (4.4%) match patient 4
    ###Check whether consistent with multiplets of other patients

    ###Find codes of multiplets that matched individual patients, later check if they match combos

    hugenvec.6multiplet.2=hugenvec26[hmatch.6multiplet.26[!is.na(hmatch.6multiplet.26)]]
    hugenvec.6multiplet.3=hugenvec36[hmatch.6multiplet.36[!is.na(hmatch.6multiplet.36)]]
    hugenvec.6multiplet.4=hugenvec46[hmatch.6multiplet.46[!is.na(hmatch.6multiplet.46)]]

    ###Now find all the combinations

    h6.combos=rep(NA,length(hugenvec26)*length(hugenvec36)+length(hugenvec26)*length(hugenvec46)+length(hugenvec36)*length(hugenvec46)+length(hugenvec26)*length(hugenvec36)*length(hugenvec46))

    count=0
    for(i in 1:length(hugenvec26))
    {
        one=as.numeric(unlist(strsplit(hugenvec26[i],split="")))
        for(j in 1:length(hugenvec36))
            {
                two=as.numeric(unlist(strsplit(hugenvec36[j],split="")))
                s=one+two
                s[s>1]=1
                count=count+1
                h6.combos[count]=paste0(s,collapse="")
                for(k in 1:length(hugenvec46))
                {
                    three=as.numeric(unlist(strsplit(hugenvec46[k],split="")))
                    s=one+three
                    s[s>1]=1
                    count=count+1
                    h6.combos[count]=paste0(s,collapse="")
                    s=two+three
                    s[s>1]=1
                    count=count+1
                    h6.combos[count]=paste0(s,collapse="")
                    s=one+two+three
                    s[s>1]=1
                    count=count+1
                    h6.combos[count]=paste0(s,collapse="")
                }
            }
    }

    h6.combos=sort(unique(h6.combos))

    ###match multiplets to combos

    hmatch.6multiplet.combos=match(genvec.6multiplet,h6.combos)

    ###346 out of 792 don't match (43.7%).  There are 91 genotypes from multiplets

    genvec.6multiplet.nomatch=genvec.6multiplet[is.na(hmatch.6multiplet.combos)]

    ###Do multiplets that match individual patients also match combos

    match.hugenvec.6multiplet.2.combos=match(hugenvec.6multiplet.2,h6.combos)
    match.hugenvec.6multiplet.3.combos=match(hugenvec.6multiplet.3,h6.combos)
    match.hugenvec.6multiplet.4.combos=match(hugenvec.6multiplet.4,h6.combos)

    ###13 of 13 multiplets that match patient 2 also match combos, 0 of 10 multiplets that match patient 3 match combos, and 0 of 35 that match patient 4 match combos

    ###Now look at imputation

    ugenvec6.unimputed=unique(dat6$genoUnimputed)
    n6.unimputed=length(ugenvec6.unimputed)
    hugenvec6.unimputed=rep(NA,n6.unimputed)
    for(i in 1:n6.unimputed) hugenvec6.unimputed[i]=paste0(unlist(strsplit(ugenvec6.unimputed[i],split=" ")),collapse="")
    which.6.9=(1:n6.unimputed)[grep(9,hugenvec6.unimputed)]
    which.6.not9=(1:n6.unimputed)[-grep(9,hugenvec6.unimputed)]
    hugenvec6.unimputed.9=hugenvec6.unimputed[which.6.9]
    hugenvec6.unimputed.not9=hugenvec6.unimputed[which.6.not9]

    ###See if multiplets match unimputed genotypes

    match.6multiplet.unimputed=match(unique(genvec.6multiplet),hugenvec6.unimputed.not9)
    count.6multiplet.unimputed=sum(is.na(match.6multiplet.unimputed))

    ###80 of 91 genotypes are not imputed, so imputation probably not a big deal


    ###
    ###
    ###SNACS7

    dat7=read.delim(paste0(dirData,"annCell_SNACS7",exptNameSuffix,"_SNACS7snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    for (k in c("genotype","genoUnimputed")) dat7[,k]=as.character(dat7[,k])
    nsnp7=length(strsplit(dat7$genotype[1],split=" ")[[1]])
    n7=nrow(dat7)
    
    #amb7=getAmbiguousGenotype(exptName="SNACS7",hashNames=unique(dat7$snacsRnd2))
    amb7=getAmbiguousGenotype(exptName="SNACS7",thresGenoCall=thresGenoCall)

    ###SNACS1, SNACS2, SNACS3, SNACS4

    dat17=read.delim(paste0(dirData,"annCell_SNACS1",exptNameSuffix,"_SNACS7snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    dat27=read.delim(paste0(dirData,"annCell_SNACS2",exptNameSuffix,"_SNACS7snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    dat37=read.delim(paste0(dirData,"annCell_SNACS3",exptNameSuffix,"_SNACS7snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    dat47=read.delim(paste0(dirData,"annCell_SNACS4",exptNameSuffix,"_SNACS7snp.txt"),sep="\t",header=TRUE,as.is=TRUE,row.names=1)
    for (k in c("genotype","genoUnimputed")) dat17[,k]=as.character(dat17[,k])
    for (k in c("genotype","genoUnimputed")) dat27[,k]=as.character(dat27[,k])
    for (k in c("genotype","genoUnimputed")) dat37[,k]=as.character(dat37[,k])
    for (k in c("genotype","genoUnimputed")) dat47[,k]=as.character(dat47[,k])
    n17=nrow(dat17)
    n27=nrow(dat27)
    n37=nrow(dat37)
    n47=nrow(dat47)

    ###genmat7 is the genotype of SNACS7

    genmat7=matrix(NA,n7,nsnp7)

    for(i in 1:n7)
    {
        genmat7[i,]=strsplit(dat7$genotype[i],split=" ")[[1]]
    }

    ###genmat 17, genmat27, genmat37 and genmat47 are genotypes of components of genmat7

    genmat17=matrix(NA,n17,nsnp7)
    for(i in 1:n17)
    {
        genmat17[i,]=strsplit(dat17$genotype[i],split=" ")[[1]]
    }

    genmat27=matrix(NA,n27,nsnp7)
    for(i in 1:n27)
    {
        genmat27[i,]=strsplit(dat27$genotype[i],split=" ")[[1]]
    }

    genmat37=matrix(NA,n37,nsnp7)
    for(i in 1:n37)
    {
        genmat37[i,]=strsplit(dat37$genotype[i],split=" ")[[1]]
    }

    genmat47=matrix(NA,n47,nsnp7)
    for(i in 1:n47)
    {
        genmat47[i,]=strsplit(dat47$genotype[i],split=" ")[[1]]
    }

    ###genvec collapses the genotypes

    genvec7=rep(NA,n7)
    for(i in 1:n7)
    {
        genvec7[i]=paste0(genmat7[i,],collapse="")
    }

    genvec17=rep(NA,n17)
    for(i in 1:n17)
    {
        genvec17[i]=paste0(genmat17[i,],collapse="")
    }

    genvec27=rep(NA,n27)
    for(i in 1:n27)
    {
        genvec27[i]=paste0(genmat27[i,],collapse="")
    }

    genvec37=rep(NA,n37)
    for(i in 1:n37)
    {
        genvec37[i]=paste0(genmat37[i,],collapse="")
    }

    genvec47=rep(NA,n47)
    for(i in 1:n47)
    {
        genvec47[i]=paste0(genmat47[i,],collapse="")
    }

    ###Make ugenvec as unique components of genvec

    ugenvec7=sort(unique(genvec7))
    ugenvec17=sort(unique(genvec17))
    ugenvec27=sort(unique(genvec27))
    ugenvec37=sort(unique(genvec37))
    ugenvec47=sort(unique(genvec47))

    ###Do the components ever match?

    match.17.27=match(ugenvec17,ugenvec27)
    match.17.37=match(ugenvec17,ugenvec37)
    match.17.47=match(ugenvec17,ugenvec47)
    match.27.37=match(ugenvec27,ugenvec37)
    match.27.47=match(ugenvec27,ugenvec47)
    match.37.47=match(ugenvec37,ugenvec47)

    ###No, components are unique

    ###Count up the calls in SNACS7

    which.71=which(dat7$snacsRnd2=="TS.1")
    which.72=which(dat7$snacsRnd2=="TS.2")
    which.73=which(dat7$snacsRnd2=="TS.3")
    which.74=which(dat7$snacsRnd2=="TS.4")
    which.712=which(dat7$snacsRnd2=="TS.1_TS.2")
    which.713=which(dat7$snacsRnd2=="TS.1_TS.3")
    which.714=which(dat7$snacsRnd2=="TS.1_TS.4")
    which.723=which(dat7$snacsRnd2=="TS.2_TS.3")
    which.724=which(dat7$snacsRnd2=="TS.2_TS.4")
    which.734=which(dat7$snacsRnd2=="TS.3_TS.4")
    which.7123=which(dat7$snacsRnd2=="TS.1_TS.2_TS.3")
    which.7124=which(dat7$snacsRnd2=="TS.1_TS.2_TS.4")
    which.7134=which(dat7$snacsRnd2=="TS.1_TS.3_TS.4")
    which.7234=which(dat7$snacsRnd2=="TS.2_TS.3_TS.4")
    which.71234=which(dat7$snacsRnd2=="TS.1_TS.2_TS.3_TS.4")

    ###4630 called 1, 560 called 2, 2719 called 3, 1415 called 4, 133 called 12, 989 called 13, 435 called 14, 4 called 23, 3 called 24, 234 called 34, 0 called 123, 0 called 124, 106 called 134, 0 called 234, 0 called 1234

    genvec.71=genvec7[which.71]
    genvec.72=genvec7[which.72]
    genvec.73=genvec7[which.73]
    genvec.74=genvec7[which.74]
    genvec.712=genvec7[which.712]
    genvec.713=genvec7[which.713]
    genvec.714=genvec7[which.714]
    genvec.723=genvec7[which.723]
    genvec.724=genvec7[which.724]
    genvec.734=genvec7[which.734]
    genvec.7123=genvec7[which.7123]
    genvec.7124=genvec7[which.7124]
    genvec.7134=genvec7[which.7134]
    genvec.7234=genvec7[which.7234]
    genvec.71234=genvec7[which.71234]
    genvec.7multiplet=genvec7[c(which.712,which.713,which.714,which.723,which.724,which.734,which.7123,which.7124,which.7134,which.7234,which.71234)]

    ###How many match expected genotypes for patient 1

    match.71.17=match(genvec.71,ugenvec17)

    ###61 of 4360 (1.4%) called patient 1 don't match

    match.71.27=match(genvec.71,ugenvec27)
    match.71.37=match(genvec.71,ugenvec37)
    match.71.47=match(genvec.71,ugenvec47)

    ###0%

    ###How many match expected genotypes for patient 2

    match.72.27=match(genvec.72,ugenvec27)

    ###5 of 560 (0.9%) called patient 2 don't match

    match.72.17=match(genvec.72,ugenvec17)
    match.72.37=match(genvec.72,ugenvec37)
    match.72.47=match(genvec.72,ugenvec47)

    ###1 out of 560 (0.2%) matches patient 1, 0 out of 560 (0%) matches patient 3, 0 out of 560 matches patient 4

    ###How many match expected genotypes for patient 3

    match.73.37=match(genvec.73,ugenvec37)

    ###44 of 2719 (1.6%) called patient 3 don't match

    match.73.17=match(genvec.73,ugenvec17)
    match.73.27=match(genvec.73,ugenvec27)
    match.73.47=match(genvec.73,ugenvec47)

    ###0%

    ###How many match expected genotypes for patient 4

    match.74.47=match(genvec.74,ugenvec47)

    ###44 of 1415 (3.1%) called patient 4 don't match

    match.74.17=match(genvec.74,ugenvec17)
    match.74.27=match(genvec.74,ugenvec27)
    match.74.37=match(genvec.74,ugenvec37)

    ###0%

    ###Now which multiplets match 2, 3 and 4

    match.7multiplet.17=match(genvec.7multiplet,ugenvec17)
    match.7multiplet.27=match(genvec.7multiplet,ugenvec27)
    match.7multiplet.37=match(genvec.7multiplet,ugenvec37)
    match.7multiplet.47=match(genvec.7multiplet,ugenvec47)

    ###48 of 1904 (2.5%) matches patient 1, 46 of 1904 (2.3%) matches patient 2, 230 of 1904 match patient 3 (12.1%), 34 of 1904 match patient 3 (1.8%)

    ###Logical ordering of genotypes

    ###Patient 1, 0EE0000E00EE00
    ###Patient 2, 0010000101E0EE
    ###Patient 3, 0010E000010100
    ###Patient 4, E00E1EE0E10100

    ###Adding an "h" for hand curated genotypes

    ###Building a function for all hand curated genotypes

    hfunction=function(string)
    {
        foo=unlist(strsplit(string,split=""))
        n=length(foo)
        which.e=which(foo=="E")
        nume=length(which.e)
        combos=2^nume
        output=rep(NA,combos)
        l01=vector("list",nume)
        reps.vector=rep(NA,nume)
        for(i in 1:nume) reps.vector[i]=2^(nume-i)
        reps.list=vector("list",nume)
        for(i in 1:nume) reps.list[[i]]=c(rep(0,reps.vector[i]),rep(1,reps.vector[i]))
        reps.mat=matrix(NA,nume,combos)
        for(i in 1:nume) reps.mat[i,]=rep(reps.list[[i]],combos/length(reps.list[[i]]))
        output.mat=matrix(NA,combos,n)
        counter=1
        for(i in 1:n)
        {
            if(foo[i]!="E") output.mat[,i]=foo[i]
            if(foo[i]=="E")
            {
                output.mat[,i]=reps.mat[counter,]
                counter=counter+1
            }
        }
        for(i in 1:combos)
        {
            output[i]=paste(output.mat[i,],collapse="")
        }
    ###    return(list(foo,n,nume,combos,output,reps.vector,reps.list,reps.mat,output.mat,output))
        output=sort(output)
        output
    }

    if (truthCallType=="manual") {
        ## NOT USED
        ## Related to truth call based on 2 rounds of SNP clustering
        ###Adding an "h" for hand curated genotypes
        hugenvec17=hfunction("0EE0000E00EE00")
        hugenvec27=hfunction("0010E000010100")
        hugenvec37=hfunction("0010000101E0EE")
        hugenvec47=hfunction("E00E1EE0E10100")
    } else {
        ## Unique genotype combinations from singlet experiments
        res=getUniqueGenoComb(exptName="SNACS7",thresGenoCall=thresGenoCall)
        hugenvec17=res[["SNACS1"]]
        hugenvec27=res[["SNACS2"]]
        hugenvec37=res[["SNACS3"]]
        hugenvec47=res[["SNACS4"]]
    }

    ###Now match to patient 1

    hmatch.71.17=match(genvec.71,hugenvec17)
    hmatch.71.27=match(genvec.71,hugenvec27)
    hmatch.71.37=match(genvec.71,hugenvec37)
    hmatch.71.47=match(genvec.71,hugenvec47)

    ###122 of 4360 (2.8%) don't match, 0 match otheres

    ###Now match to patient 2

    hmatch.72.17=match(genvec.72,hugenvec17)
    hmatch.72.27=match(genvec.72,hugenvec27)
    hmatch.72.37=match(genvec.72,hugenvec37)
    hmatch.72.47=match(genvec.72,hugenvec47)

    ###17 of 560 (3.0%) don't match, then 1 0 1

    ###Now match to patient 3

    hmatch.73.17=match(genvec.73,hugenvec17)
    hmatch.73.27=match(genvec.73,hugenvec27)
    hmatch.73.37=match(genvec.73,hugenvec37)
    hmatch.73.47=match(genvec.73,hugenvec47)

    ###91 of 2719 (3.3%) don't match, then 0 0 0

    ###Now match to patient 4

    hmatch.74.17=match(genvec.74,hugenvec17)
    hmatch.74.27=match(genvec.74,hugenvec27)
    hmatch.74.37=match(genvec.74,hugenvec37)
    hmatch.74.47=match(genvec.74,hugenvec47)

    ###48 of 1415 (3.4%) don't match, then 0, 0, 0

    ###Now which multiplets match 1, 2, 3 and 4

    hmatch.7multiplet.17=match(genvec.7multiplet,hugenvec17)
    hmatch.7multiplet.27=match(genvec.7multiplet,hugenvec27)
    hmatch.7multiplet.37=match(genvec.7multiplet,hugenvec37)
    hmatch.7multiplet.47=match(genvec.7multiplet,hugenvec47)

    ###11 (0.6%) match patient1,  3 (0.2%) match patient 2, 60 (3.2%) match patient 3, and 20 (1.1%) match patient 4
    ###Check whether consistent with multiplets of other patients

    ###Find codes of multiplets that matched individual patients, later check if they match combos

    hugenvec.7multiplet.1=hugenvec17[hmatch.7multiplet.17[!is.na(hmatch.7multiplet.17)]]
    hugenvec.7multiplet.2=hugenvec27[hmatch.7multiplet.27[!is.na(hmatch.7multiplet.27)]]
    hugenvec.7multiplet.3=hugenvec37[hmatch.7multiplet.37[!is.na(hmatch.7multiplet.37)]]
    hugenvec.7multiplet.4=hugenvec47[hmatch.7multiplet.47[!is.na(hmatch.7multiplet.47)]]

    ###Now find all the combinations

    h7.combos=NULL

    for(i in 1:length(hugenvec17))
    {
        one=as.numeric(unlist(strsplit(hugenvec17[i],split="")))
        for(j in 1:length(hugenvec27))
        {
            two=as.numeric(unlist(strsplit(hugenvec27[j],split="")))
            s=one+two
            s[s>1]=1
            h7.combos=c(h7.combos,paste0(s,collapse=""))
            for(k in 1:length(hugenvec37))
            {
                three=as.numeric(unlist(strsplit(hugenvec37[k],split="")))
                s=one+three
                s[s>1]=1
                h7.combos=c(h7.combos,paste0(s,collapse=""))
                s=two+three
                s[s>1]=1
                h7.combos=c(h7.combos,paste0(s,collapse=""))
                s=one+two+three
                s[s>1]=1
                h7.combos=c(h7.combos,paste0(s,collapse=""))
                for(l in 1:length(hugenvec47))
                {
                    four=as.numeric(unlist(strsplit(hugenvec47[l],split="")))
                    s=one+four
                    s[s>1]=1
                    h7.combos=c(h7.combos,paste0(s,collapse=""))
                    s=two+four
                    s[s>1]=1
                    h7.combos=c(h7.combos,paste0(s,collapse=""))
                    s=three+four
                    s[s>1]=1
                    h7.combos=c(h7.combos,paste0(s,collapse=""))
                    s=one+two+four
                    s[s>1]=1
                    h7.combos=c(h7.combos,paste0(s,collapse=""))
                    s=one+three+four
                    s[s>1]=1
                    h7.combos=c(h7.combos,paste0(s,collapse=""))
                    s=two+three+four
                    s[s>1]=1
                    h7.combos=c(h7.combos,paste0(s,collapse=""))
                    s=one+two+three+four
                    s[s>1]=1
                    h7.combos=c(h7.combos,paste0(s,collapse=""))
                }
            }
        }
    }

    h7.combos=sort(unique(h7.combos))

    ###match multiplets to combos

    hmatch.7multiplet.combos=match(genvec.7multiplet,h7.combos)

    ###903 out of 1904 don't match (47.4%).  There are 576 genotypes from multiplets

    genvec.7multiplet.nomatch=genvec.7multiplet[is.na(hmatch.7multiplet.combos)]

    ###Do multiplets that match individual patients also match combos

    match.hugenvec.7multiplet.1.combos=match(hugenvec.7multiplet.1,h7.combos)
    match.hugenvec.7multiplet.2.combos=match(hugenvec.7multiplet.2,h7.combos)
    match.hugenvec.7multiplet.3.combos=match(hugenvec.7multiplet.3,h7.combos)
    match.hugenvec.7multiplet.4.combos=match(hugenvec.7multiplet.4,h7.combos)

    ###0 of 11 multiplets that match patient 1 also match combos, 3 of 13 multiplets that match patient 2 match combos, and 60 of 60 that match patient 3 match combos, 4 of 20 multiplets that match patient 4 match combos


    combos.2=function(single.list)
    {
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

    combos.3=function(single.list)
    {
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

    combos.4=function(single.list)
    {
        one=single.list[[1]]
        two=single.list[[2]]
        three=single.list[[3]]
        four=single.list[[4]]
        output=NULL

        for(i in 1:length(one))
        {
            seq1=as.numeric(unlist(strsplit(one[i],split="")))
            for(j in 1:length(two))
            {
                seq2=as.numeric(unlist(strsplit(two[j],split="")))
                s=seq1+seq2
                s[s>1]=1
                output=c(output,paste0(s,collapse=""))
                for(k in 1:length(three))
                {
                    seq3=as.numeric(unlist(strsplit(three[k],split="")))
                    s=seq1+seq3
                    s[s>1]=1
                    output=c(output,paste0(s,collapse=""))
                    s=seq2+seq3
                    s[s>1]=1
                    output=c(output,paste0(s,collapse=""))
                    s=seq1+seq2+seq3
                    s[s>1]=1
                    output=c(output,paste0(s,collapse=""))
                    for(l in 1:length(four))
                    {
                        seq4=as.numeric(unlist(strsplit(four[l],split="")))
                        s=seq1+seq4
                        s[s>1]=1
                        output=c(output,paste0(s,collapse=""))
                        s=seq2+seq4
                        s[s>1]=1
                        output=c(output,paste0(s,collapse=""))
                        s=seq3+seq4
                        s[s>1]=1
                        output=c(output,paste0(s,collapse=""))
                        s=seq1+seq2+seq4
                        s[s>1]=1
                        output=c(output,paste0(s,collapse=""))
                        s=seq1+seq3+seq4
                        s[s>1]=1
                        output=c(output,paste0(s,collapse=""))
                        s=seq2+seq3+seq4
                        s[s>1]=1
                        output=c(output,paste0(s,collapse=""))
                        s=seq1+seq2+seq3+seq4
                        s[s>1]=1
                        output=c(output,paste0(s,collapse=""))
                    }
                }
            }
        }
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

    truth.function=function(gv,single.list,amb=NULL)
    {
        #gv=genvec6; single.list=list(hugenvec26,hugenvec36,hugenvec46); amb=amb6
        #gv=genvec7; single.list=list(hugenvec17,hugenvec27,hugenvec37,hugenvec47); amb=amb7
        
        cat("\n------------ truth.function: SNACS",length(single.list)+3," --------------\n",sep="")
        timeStamp=Sys.time()
        
    ###How many samples have been multiplexed
        n.samples=length(single.list)
        if(n.samples==2) combos=combos.2(single.list)
        if(n.samples==3) combos=combos.3(single.list)
        if(n.samples==4) combos=combos.4(single.list)
    ###    if(n.samples==3) combos=combos.3(single.list)
    ###    if(n.samples==4) combos=combos.4(single.list)
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
        combos.matrix=matrix(NA,length(combos),length(unlist(strsplit(combos[1],split=""))))
        for(i in 1:n.combos) combos.matrix[i,]=unlist(strsplit(combos[i],split=""))
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
        output=list(ugv=ugv,match.singles.exact=match.singles.exact,match.combos.exact=match.combos.exact,no.matches=no.matches,singles=singles,singles.matrix=singles.matrix,combos.matrix=combos.matrix,match.combos.extended=match.combos.extended,indicator.multiplet=indicator.multiplet,match.gv.ugv=match.gv.ugv,matches.vector=matches.vector,final.calls=final.calls)
        
        timeStamp=c(timeStamp,Sys.time())
        print(timeStamp)
        print(diff(timeStamp))

        output
        final.calls
    }
    
    foo=truth.function(genvec5,list(hugenvec15,hugenvec25),amb=amb5)
    
    #makeHashPairScatterPlot(exptName="SNACS5",exptNameSuffix)
    if (!file.exists(dirOutput)) dir.create(file.path(dirOutput))
    tbl=data.frame(id=rownames(dat5))
    tbl$hashCall_truth=foo
    write.table(tbl,file=paste0(dirOutput,"trueHashCall_SNACS5",exptNameSuffix,"_SNACS5snp.txt"),col.names=T,row.names=F, sep="\t",quote=F)

    foo2=truth.function(genvec6,list(hugenvec26,hugenvec36,hugenvec46),amb=amb6)
    
    #makeHashPairScatterPlot(exptName="SNACS6",exptNameSuffix)
    tbl=data.frame(id=rownames(dat6))
    tbl$hashCall_truth=foo2
    write.table(tbl,file=paste0(dirOutput,"trueHashCall_SNACS6",exptNameSuffix,"_SNACS6snp.txt"),col.names=T,row.names=F, sep="\t",quote=F)

    foo3=truth.function(genvec7,list(hugenvec17,hugenvec27,hugenvec37,hugenvec47),amb=amb7)

    #makeHashPairScatterPlot(exptName="SNACS7",exptNameSuffix)
    tbl=data.frame(id=rownames(dat7))
    tbl$hashCall_truth=foo3
    write.table(tbl,file=paste0(dirOutput,"trueHashCall_SNACS7",exptNameSuffix,"_SNACS7snp.txt"),col.names=T,row.names=F, sep="\t",quote=F)

    timeStamp=Sys.time()
    print(timeStamp)

}

####################################################################
####################################################################
