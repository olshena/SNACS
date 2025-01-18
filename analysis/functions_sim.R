
## ------------------------------
modifyCallName=function(snacsCall) {
    snacsCallThis=gsub("TS.","",snacsCall,fixed=T)
    snacsCallThis=gsub(",","_",snacsCallThis,fixed=T)
    ifelse(length(grep("_",snacsCallThis))!=0 | snacsCallThis=="Doublet","multiplet",paste0("TS.",snacsCallThis))
}


truthName2snacsCallName=function(truth,hashNameInit=NULL,hashNameFinal=NULL) {
    if (is.null(hashNameInit)) {
        truthThis=gsub(",","_TS.",paste0("TS.",truth),fixed=T)
    } else {
        truthThis=sapply(truth,function(truth,hashNameInit,hashNameFinal) {
                truthThis=strsplit(truth,",")[[1]]
                j=match(truthThis,hashNameInit); j1=which(!is.na(j)); j2=j[j1]
                truthThis=paste(hashNameFinal[j2],collapse="_")
                truthThis
            },hashNameInit,hashNameFinal,USE.NAMES=F)
    }
    truthThis
}

getColVarInfo=function(annCell,annHash) {
    cellAnnoVar=c("truth","snacsPlusDoubletD","doubletD","snacsRnd2","snacsRnd1","clustBestSNPs_hclust")
    cellAnnoName=c("truth","snacs+DD","doubletD","snacsRnd2","snacsRnd1","bestSNPsCluster")
    cellAnnoVar=c("truth","snacsPlusDoubletD","doubletD","snacsRnd2","clustSplitRnd2","snacsRnd1","clustSplitRnd1","clustBestSNPs_hclust")
    cellAnnoName=c("truth","snacs+DD","doubletD","snacsRnd2","clustSplitRnd2","snacsRnd1","clustSplitRnd1","bestSNPsCluster")
    k=which(cellAnnoVar%in%names(annCell))
    cellAnnoVar=cellAnnoVar[k]; cellAnnoName=cellAnnoName[k]
    cellAnnoVar=c(cellAnnoVar,annHash$hashNames)
    cellAnnoName=c(cellAnnoName,annHash$hashNames)
    col_var_info=list()
    if ("truth"%in%names(annCell)) {
        multNames=unique(annCell$truth)
        col_var_info[["truth"]]=list(level=unique(annCell$truth))
        multNames=sort(unique(multNames[!multNames%in%annHash$hashNames]))
        colVec=annHash$hashColors
        if (length(multNames)!=0) colVec=c(colVec,grDevices::gray((1:length(multNames))/length(multNames)))
        multNames=c(annHash$hashNames,multNames)
        for (vId in which(names(col_var_info)%in%intersect(c("truth"),names(annCell)))) {
            k=match(col_var_info[[vId]]$level,multNames)
            col_var_info[[vId]]$level=multNames[k]
            col_var_info[[vId]]$color=colVec[k]
        }
    }
    invisible(list(cellAnnoVar=cellAnnoVar,cellAnnoName=cellAnnoName,col_var_info=col_var_info))
}


createHeatmapForSim=function(snacsObj,colVarInfo=NULL,cell_anno_var,cell_anno_name=NULL,col_dend=T,row_dend=F,h_title=NULL,outputFileName="heatmap",outputFormat=c("","pdf","png","none"),write2Table=F) {
    #colVarInfo=NULL; cell_anno_var; cell_anno_name=NULL; col_dend=T; row_dend=F; h_title=NULL; outputFileName="heatmap"; outputFormat=c("","pdf","png","none"); write2Table=F
    #colVarInfo=col_var_info; cell_anno_var=cellAnnoVar; cell_anno_name=cellAnnoName; col_dend=TRUE; row_dend=FALSE; outputFileName=paste0("heatmap_",snacsObj$exptName); outputFormat="pdf"
    
    ################################################
    ## Parameters

    outputFormat=outputFormat[1]
    
    snp_anno_var=snp_anno_name=NULL

    ################################################
    ## Parameters - optional

    ## Whether to generate sample & heatmap color legends
    showLegend=F

    ################################################
    ## Compute cosine cross-distances between the rows of matrices
    getSKmeansDist=function(x) {stats::as.dist(skmeans::skmeans_xdist(x))}
    distfun=getSKmeansDist
    linkMethod="ward.D2"

    ################################################
    if (is.null(snp_anno_var)) hasRowAnno=F else hasRowAnno=T
    if (is.null(cell_anno_name)) cell_anno_name=cell_anno_var
    if (is.null(snp_anno_name)) snp_anno_name=snp_anno_var
    cell_anno_name=paste0(cell_anno_name," ") ## Cell color bar names
    snp_anno_name=paste0(snp_anno_name," ") ## SNP color bar names
    k=match(c("SnpEff_Annotation_Impact","meanQuality","meanDepth"),snp_anno_var); k1=which(!is.na(k)); k2=k[k1]
    if (length(k1)!=0) snp_anno_name[k2]=paste0(c("Impact","Quality","Depth")[k1]," ")

    snacsObj$annCell=cbind(snacsObj$annCell[,which(!names(snacsObj$annCell)%in%rownames(snacsObj$hashes))],t(snacsObj$hashes),stringsAsFactors=F)
    
    fName=outputFileName

    ################################################
    ## Folder where heatmap will be generated
    
    if (!outputFormat%in%c("","none")) {
    #if (outputFormat!="") {
        subDir="../output/"; if (!file.exists(subDir)) dir.create(file.path(subDir))
        subDir="../output/heatmap/"; if (!file.exists(subDir)) dir.create(file.path(subDir))
    }

    ################################################
    ## "generate_heatmap()"  parameters
    
    if (col_dend) ncc=nrow(snacsObj$annHash) else ncc=NA; if (!row_dend) ncr=NA else ncr=NA
    heatmap_color=c("navy","gray85","navy")
    limHM=c(0,1)
    plotInfo=list(margins=c(5,.2),cexRow=0.2,cexCol=0.2,cexRowSide=1,cexColSide=.7,colorCatCol=c("skyblue", "blue", "yellow", "purple", "black", "red", "orange", "green", "cyan", "darkgreen","grey","brown","pink","salmon","palegreen"))
    if (outputFormat=="png") plotInfo$cexColSide=1.5

    if (is.null(colVarInfo)) {col_var_info=list()} else {col_var_info=colVarInfo}
    
    ## Set colors for the SNACS calls
    ## Multiplets are shades of gray
    multNames=""
    if (T) {
        ## devtools::check("SNACS") gives error with grDevices::gray() if section is out of the if () statement
        if ("snacsRnd1"%in%names(snacsObj$annCell)) {
            multNames=unique(snacsObj$annCell$snacsRnd1)
            col_var_info[["snacsRnd1"]]=list(level=unique(snacsObj$annCell$snacsRnd1))
        }
        if ("snacsRnd2"%in%names(snacsObj$annCell)) {
            multNames=c(multNames,unique(snacsObj$annCell$snacsRnd2))
            col_var_info[["snacsRnd2"]]=list(level=unique(snacsObj$annCell$snacsRnd2))
        }
        if ("snacsPlusDoubletD"%in%names(snacsObj$annCell)) {
            multNames=c(multNames,unique(snacsObj$annCell$snacsPlusDoubletD))
            col_var_info[["snacsPlusDoubletD"]]=list(level=unique(snacsObj$annCell$snacsPlusDoubletD))
        }
        multNames=sort(unique(multNames[!multNames%in%snacsObj$annHash$hashNames]))
    }
    colVec=snacsObj$annHash$hashColors
    if (length(multNames)!=0) colVec=c(colVec,grDevices::gray((1:length(multNames))/length(multNames)))
    multNames=c(snacsObj$annHash$hashNames,multNames)
    if ("doubletD"%in%names(snacsObj$annCell)) {
        col_var_info[["doubletD"]]=list(color=c("navy","white"),level=c("Singlet","Doublet"))
    }

    row_var_info=list(SnpEff_Annotation_Impact=list(color=c("brown","yellow","orange","red","white"),
                            level=c("HIGH","LOW","MODERATE","MODIFIER","")),
                            meanDepth=list(limit=c(20,80)),meanQuality=list(limit=c(25,95)))

    for (vId in which(names(col_var_info)%in%intersect(c("snacsRnd1","snacsRnd2","snacsPlusDoubletD"),names(snacsObj$annCell)))) {
        k=match(col_var_info[[vId]]$level,multNames)
        col_var_info[[vId]]$level=multNames[k]
        col_var_info[[vId]]$color=colVec[k]
    }

    for (vId in which(names(col_var_info)%in%names(snacsObj$annCell))) {
        k=match(sort(unique(snacsObj$annCell[,names(col_var_info)[vId]])),col_var_info[[vId]]$level)
        col_var_info[[vId]]$level=col_var_info[[vId]]$level[k]
        col_var_info[[vId]]$color=col_var_info[[vId]]$color[k]
    }
    for (vId in which(names(row_var_info)%in%names(snacsObj$annSNP))) {
        k=match(sort(unique(snacsObj$annSNP[,names(row_var_info)[vId]])),row_var_info[[vId]]$level)
        row_var_info[[vId]]$level=row_var_info[[vId]]$level[k]
        row_var_info[[vId]]$color=row_var_info[[vId]]$color[k]
    }

    ################################################
    if (!outputFormat%in%c("","none")) {
        dirH="../output/heatmap/"
    }

    switch(outputFormat,
        "png"={grDevices::png(paste(dirH,fName,"_%1d.png",sep=""),width=480*2,height=480*2)},
        "pdf"={grDevices::pdf(paste(dirH,fName,".pdf",sep=""),width=7,height=7)}
    )
    
    ## Generate heatmap
    if (outputFormat=="none") {
        distMat=distfun(t(snacsObj$mut))
        clustObj=list(colClust=stats::hclust(distMat,method=linkMethod))
    } else {
        clustObj=heatmap4::generate_heatmap(x=snacsObj$mut,distfun=distfun,methodR=linkMethod,methodC=linkMethod,col_anno=T,row_anno=hasRowAnno,col_info=snacsObj$annCell,row_info=snacsObj$annSNP,col_dend=col_dend,row_dend=row_dend,col_lab=F,row_lab=F,zlm=limHM,ncc=ncc,ncr=ncr,heatmap_color=heatmap_color,col_anno_var=cell_anno_var,col_anno_name=cell_anno_name,col_var_info=col_var_info,row_anno_var=snp_anno_var,row_anno_name=snp_anno_name,row_var_info=row_var_info,densColor=NULL,h_title=h_title,plot_info=plotInfo,input_legend=showLegend)
    }
    
    if (!outputFormat%in%c("","none")) grDevices::dev.off()

    if (write2Table) {
        ## Write the cell & SNP annotation tables ordered as in heatmap
        tbl=heatmap4::cutCluster(clustObj$rowClust,ann=snacsObj$annSNP,nClust=ncr,rev=T)
        utils::write.table(tbl,paste0(dirH,"clustInfoSnp_",fName,".txt"),sep="\t",col.names=T,row.names=F,quote=F)
        tbl=heatmap4::cutCluster(clustObj$colClust,ann=snacsObj$annCell,nClust=ncc,rev=F)
        utils::write.table(tbl,paste0(dirH,"clustInfoCell_",fName,".txt"),sep="\t",col.names=T,row.names=F,quote=F)
    }
    
    invisible(clustObj)
}

plotHashPair=function(so,callCol="snacsRnd2",dirOutput="../output/") {
    dirOutput2=paste0(dirOutput,"plotHashPair/")
    if (!file.exists(dirOutput2)) dir.create(file.path(dirOutput2))

    lim=range(c(so$hashes),na.rm=T)

    #png(paste0(dirOutput2,"plot_hashPair_",so$exptName,".png"),width=2*270,height=2*270)

    png(paste0(dirOutput2,"plot_hashPair_",so$exptName,".png"),width=3*270,height=2*270)
    par(mfrow=c(2,3))

    par(mar=c(5, 5, 5, 3) + 0.1)
    
    #if (callCol=="truth") {callThis=paste0("sample",so$annCell$truth)} else {callThis=so$annCell[,callCol]}
    callThis=sapply(so$annCell[,callCol],modifyCallName,USE.NAMES=F)

    colVec=rep("black",nrow(so$annCell))
    j=match(callThis,so$annHash$hashNames); j1=which(!is.na(j)); j2=j[j1]
    colVec[j1]=so$annHash$hashColors[j2]
    
    for (hId1 in 1:(nrow(so$annHash)-1)) {
        if (callCol=="truth") {
            #j1=which(so$annCell$truth==hId1)
            j1=which(so$annCell[,callCol]==so$annHash$hashNames[hId1])
        } else {
            j1=which(so$annCell[,callCol]==so$annHash$hashNames[hId1])
        }
        for (hId2 in (hId1+1):nrow(so$annHash)) {
            if (callCol=="truth") {
                #j2=which(so$annCell$truth==hId2)
                j2=which(so$annCell[,callCol]==so$annHash$hashNames[hId2])
            } else {
                j2=which(so$annCell[,callCol]==so$annHash$hashNames[hId2])
            }

            if (F) {
                colVec=rep("grey",nrow(so$annCell))
                colVec[j1]=so$annHash$hashColors[hId1]
                colVec[j2]=so$annHash$hashColors[hId2]
                colVec[-c(j1,j2)]="black"
            }

            hashNames=sub("TS.","Sample",so$annHash$hashNames,fixed=T)
            header=so$exptName
            header=sub("samples_","samples\n",so$exptName)
            plot(so$hashes[hId1,],so$hashes[hId2,],xlim=lim,ylim=lim,main=header,xlab=paste0(hashNames[hId1],": Hash"),ylab=paste0(hashNames[hId2],": Hash"),cex.main=2,cex.lab=2,cex.axis=1.5,col=colVec,pch=16)
            j=which(colVec=="black")
            if (length(j)!=0) {
                points(so$hashes[hId1,j],so$hashes[hId2,j],pch=16)
            }
        }
    }
    dev.off()
}
