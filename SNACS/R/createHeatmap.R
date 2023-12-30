#' Generate heatmap of mutation data
#'
#' @param snacsObj SNACSList object
#' @param cell_anno_var Character vector. Cell variables for which color bars are to be added to the heatmap. Default is NULL
#' @param cell_anno_name  Character vector. Names of column color bars. Default is NULL
#' @param col_dend Logical. Displays cell dendrogram. Default is TRUE
#' @param row_dend Logical. Displays SNP dendrogram. Default is FALSE
#' @param h_title Character. Gives the heatmap output a title. Default is NULL
#' @param outputFileName Character. Output file name. Default is "heatmap"
#' @param outputFormat Character. Output file type. Default is "" which outputs to the standard output
#' @param write2Table Logical. Write the cell & SNP annotation tables ordered as in heatmap. Default is FALSE

#' @return A SNACSList object
#' @export
createHeatmap=function(snacsObj,cell_anno_var,cell_anno_name=NULL,col_dend=T,row_dend=F,h_title=NULL,outputFileName="heatmap",outputFormat=c("","pdf","png"),write2Table=F) {
    
    ################################################
    ## Parameters

    outputFormat=outputFormat[1]
    
    snp_anno_var=snp_anno_name=NULL

    ################################################
    ## Parameters - optional

    ## Whether to generate the heatmap in a subfolder
    subsetCellFlag=""

    ## Whether to generate sample & heatmap color legends
    showLegend=T
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
    
    if (outputFormat!="") {
        subDir="../output/"; if (!file.exists(subDir)) dir.create(file.path(subDir))
        subDir="../output/heatmap/"; if (!file.exists(subDir)) dir.create(file.path(subDir))
    }

    ################################################
    ## "generate_heatmap()"  parameters
    
    if (col_dend) ncc=nrow(snacsObj$annHash) else ncc=NA; if (!row_dend) ncr=NA else ncr=NA
    #heatmap_color=c("orangered1","royalblue3","white")
    #heatmap_color=c("navy","gray85","orangered1")
    heatmap_color=c("navy","gray85","navy")
    limHM=c(0,1)
    plotInfo=list(margins=c(5,.2),cexRow=0.2,cexCol=0.2,cexRowSide=1,cexColSide=.7,colorCatCol=c("skyblue", "blue", "yellow", "purple", "black", "red", "orange", "green", "cyan", "darkgreen","grey","brown","pink","salmon","palegreen"))
    if (outputFormat=="png") plotInfo$cexColSide=1.5

    col_var_info=list()
    
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
    colVec=c(snacsObj$annHash$hashColors,grDevices::gray((1:length(multNames))/length(multNames)))
    multNames=c(snacsObj$annHash$hashNames,multNames)
    if ("doubletD"%in%names(snacsObj$annCell)) {
        col_var_info[["doubletD"]]=list(color=c("navy","white"),level=c("Singlet","Doublet"))
        #k=which(col_var_info[["doubletD"]]$level%in%snacsObj$annCell$hashCall_HTOdemux)
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
    if (outputFormat!="") {
        subDir="../output/heatmap/"; if (!file.exists(subDir)) dir.create(file.path(subDir))
        if (subsetCellFlag!="") {subDir=paste0(subDir,subsetCellFlag,"/"); if (!file.exists(subDir)) dir.create(file.path(subDir))}
        dirH=subDir
    }

    switch(outputFormat,
        "png"={grDevices::png(paste(dirH,fName,"_%1d.png",sep=""),width=480*2,height=480*2)},
        "pdf"={grDevices::pdf(paste(dirH,fName,".pdf",sep=""),width=7,height=7)}
    )
    
    ## Generate heatmap
    clustObj=heatmap4::generate_heatmap(x=snacsObj$mut,distfun=distfun,methodR=linkMethod,methodC=linkMethod,col_anno=T,row_anno=hasRowAnno,col_info=snacsObj$annCell,row_info=snacsObj$annSNP,col_dend=col_dend,row_dend=row_dend,col_lab=F,row_lab=F,zlm=limHM,ncc=ncc,ncr=ncr,heatmap_color=heatmap_color,col_anno_var=cell_anno_var,col_anno_name=cell_anno_name,col_var_info=col_var_info,row_anno_var=snp_anno_var,row_anno_name=snp_anno_name,row_var_info=row_var_info,densColor=NULL,h_title=h_title,plot_info=plotInfo,input_legend=showLegend)
    
    if (outputFormat!="") grDevices::dev.off()
    
    if (write2Table) {
        ## Write the cell & SNP annotation tables ordered as in heatmap
        tbl=heatmap4::cutCluster(clustObj$rowClust,ann=snacsObj$annSNP,nClust=ncr,rev=T)
        utils::write.table(tbl,paste0(dirH,"clustInfoSnp_",fName,".txt"),sep="\t",col.names=T,row.names=F,quote=F)
        tbl=heatmap4::cutCluster(clustObj$colClust,ann=snacsObj$annCell,nClust=ncc,rev=F)
        utils::write.table(tbl,paste0(dirH,"clustInfoCell_",fName,".txt"),sep="\t",col.names=T,row.names=F,quote=F)
    }
    
    invisible(clustObj)
}

###########################################################
###########################################################
