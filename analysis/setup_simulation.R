## Create initial data from single-sample experiments 1-4 that will be used to simulate multi-sample experiments

####################################################################
####################################################################

library(SNACS)

## ------------------------------
## Extract raw data from hdf5 file. Then create SNACSList objects

hashColors <- c("indianred2","green3","dodgerblue3","magenta3")
dirDataRaw <- "../data/data_sequenced/"
dirData <- "../data/"

for (exptName in paste0("SNACS",1:7)) {
    cat("\n\n------------- ",exptName,"\n")
    switch(exptName,
        "SNACS1"={
            fileName="GSM8066757.hdf5"
            hashNames=c("TS.1","TS.2")
        },
        "SNACS2"={
            fileName="GSM8066758.hdf5"
            hashNames=c("TS.2","TS.3")
        },
        "SNACS3"={
            fileName="GSM8066759.hdf5"
            hashNames=c("TS.3","TS.4")
        },
        "SNACS4"={
            fileName="GSM8066760.hdf5"
            hashNames=c("TS.1","TS.4")
        },
        "SNACS5"={
            fileName="GSM8066761.hdf5"
            hashNames=c("TS.1","TS.2")
        },
        "SNACS6"={
            fileName="GSM8066762.hdf5"
            hashNames=c("TS.2","TS.3","TS.4")
        },
        "SNACS7"={
            fileName="GSM8066763.hdf5"
            hashNames=c("TS.1","TS.2","TS.3","TS.4")
        }
    )
    h5toList=h5readForSNACS(file=paste0(dirDataRaw,fileName))
    
    snacsObj=SNACSList(mut=h5toList$mut,hashes=h5toList$hashes[hashNames,],exptName=exptName,hashColors=hashColors[match(hashNames,c("TS.1","TS.2","TS.3","TS.4"))],
                          depthTotal=h5toList$depthTotal,depthAlt=h5toList$depthAlt,annCell=h5toList$annCell,annSNP=h5toList$annSNP)
    rm(h5toList)
    snacsObj$annSNP$desc=getSNPdesc(snacsObj$annSNP$desc)
    save(snacsObj,file=paste0(dirData,"snacsObj_init_",exptName,"_unfilt.RData"))
}

## ------------------------------
###Load the four individual samples

load("../data/snacsObj_init_SNACS1_unfilt.RData")
s1=snacsObj

load("../data/snacsObj_init_SNACS2_unfilt.RData")
s2=snacsObj

load("../data/snacsObj_init_SNACS3_unfilt.RData")
s3=snacsObj

load("../data/snacsObj_init_SNACS4_unfilt.RData")
s4=snacsObj

rm(snacsObj)

## ------------------------------
###mut: raw mutation data on 0,1 scale, #variant by #cells
###hashes: hash value, 2 by #cells, each row a different hash
###exptName: experiments name, like SNACS1
###depthTotal: depth total count, min is 0, #variant by #cells
###depthAlt: depth total alternative, min is -1, check on this
###annHash: The two hashes used, like TS.1 is indianred2, TS.2 is green3
###annCell: cell barcode, #cells by 2 (cellnumber,barcode)
###annSNP: variant annotations, #variants by 2 (snpnumber, snp annotation)
###processLevels: process level, in this case "raw"

####Number of variants

nv1=nrow(s1$mut)
nv2=nrow(s2$mut)
nv3=nrow(s3$mut)
nv4=nrow(s4$mut)

####Number of cells

nc1=ncol(s1$mut)
nc2=ncol(s2$mut)
nc3=ncol(s3$mut)
nc4=ncol(s4$mut)

###Find the set of common variants and make a reduced keeper data set from it

match.s1.s2.variants=match(s1$annSNP[,2],s2$annSNP[,2])
match.s1.s3.variants=match(s1$annSNP[,2],s3$annSNP[,2])
match.s1.s4.variants=match(s1$annSNP[,2],s4$annSNP[,2])

variantkeepers=which(!is.na(match.s1.s2.variants) & !is.na(match.s1.s3.variants) & !is.na(match.s1.s4.variants))
variantkeepers1=variantkeepers
variantkeepers2=match.s1.s2.variants[variantkeepers]
variantkeepers3=match.s1.s3.variants[variantkeepers]
variantkeepers4=match.s1.s4.variants[variantkeepers]

annSNP=s1$annSNP[variantkeepers1,]
annSNP$id=paste0("snp",1:nrow(annSNP))
save(annSNP,file="annSNP_forSimulation.RData")

s1$variantkeepers=variantkeepers1
s2$variantkeepers=variantkeepers2
s3$variantkeepers=variantkeepers3
s4$variantkeepers=variantkeepers4

s1$mutkeepers=s1$mut[s1$variantkeepers,]
s2$mutkeepers=s2$mut[s2$variantkeepers,]
s3$mutkeepers=s3$mut[s3$variantkeepers,]
s4$mutkeepers=s4$mut[s4$variantkeepers,]

s1$annsnpkeepers=s1$mut[s1$variantkeepers,]
s2$annsnpkeepers=s2$mut[s2$variantkeepers,]
s3$annsnpkeepers=s3$mut[s3$variantkeepers,]
s4$annsnpkeepers=s4$mut[s4$variantkeepers,]

s1$depthtotalvariantkeepers=s1$depthTotal[s1$variantkeepers,]
s2$depthtotalvariantkeepers=s2$depthTotal[s2$variantkeepers,]
s3$depthtotalvariantkeepers=s3$depthTotal[s3$variantkeepers,]
s4$depthtotalvariantkeepers=s4$depthTotal[s4$variantkeepers,]

s1$depthaltvariantkeepers=s1$depthAlt[s1$variantkeepers,]
s2$depthaltvariantkeepers=s2$depthAlt[s2$variantkeepers,]
s3$depthaltvariantkeepers=s3$depthAlt[s3$variantkeepers,]
s4$depthaltvariantkeepers=s4$depthAlt[s4$variantkeepers,]

nkeepers=length(variantkeepers)

###Indices of higher hash
###For S1, S2 and S3 the proper hash is when higher is TRUE, for S4 it is higher is FALSE

higher1=rep(TRUE,nc1)
higher2=rep(TRUE,nc2)
higher3=rep(TRUE,nc3)
higher4=rep(TRUE,nc4)

higher1[which(s1$hashes[2,]>s1$hashes[1,])]=FALSE
higher2[which(s2$hashes[2,]>s2$hashes[1,])]=FALSE
higher3[which(s3$hashes[2,]>s3$hashes[1,])]=FALSE
higher4[which(s4$hashes[2,]>s4$hashes[1,])]=FALSE

s1$hashkeepers=which(higher1)
s2$hashkeepers=which(higher2)
s3$hashkeepers=which(higher3)
s4$hashkeepers=which(!higher4)

s1$mutkeepershashkeepers=s1$mutkeepers[,s1$hashkeepers]
s2$mutkeepershashkeepers=s2$mutkeepers[,s2$hashkeepers]
s3$mutkeepershashkeepers=s3$mutkeepers[,s3$hashkeepers]
s4$mutkeepershashkeepers=s4$mutkeepers[,s4$hashkeepers]

s1$hasheskeepers=s1$hashes[,s1$hashkeepers]
s2$hasheskeepers=s2$hashes[,s2$hashkeepers]
s3$hasheskeepers=s3$hashes[,s3$hashkeepers]
s4$hasheskeepers=s4$hashes[,s4$hashkeepers]

s1$anncellkeepers=s1$annCell[s1$hashkeepers,]
s2$anncellkeepers=s2$annCell[s2$hashkeepers,]
s3$anncellkeepers=s3$annCell[s3$hashkeepers,]
s4$anncellkeepers=s4$annCell[s4$hashkeepers,]

s1$depthtotalvariantkeepershashkeepers=s1$depthtotalvariantkeepers[,s1$hashkeepers]
s2$depthtotalvariantkeepershashkeepers=s2$depthtotalvariantkeepers[,s2$hashkeepers]
s3$depthtotalvariantkeepershashkeepers=s3$depthtotalvariantkeepers[,s3$hashkeepers]
s4$depthtotalvariantkeepershashkeepers=s4$depthtotalvariantkeepers[,s4$hashkeepers]

s1$depthaltvariantkeepershashkeepers=s1$depthaltvariantkeepers[,s1$hashkeepers]
s2$depthaltvariantkeepershashkeepers=s2$depthaltvariantkeepers[,s2$hashkeepers]
s3$depthaltvariantkeepershashkeepers=s3$depthaltvariantkeepers[,s3$hashkeepers]
s4$depthaltvariantkeepershashkeepers=s4$depthaltvariantkeepers[,s4$hashkeepers]

###Harmonize the means and variances
###Proper is hash for labeled cell, improper is not that hash

mean.s1hashkeepers=apply(s1$hasheskeepers,1,mean)
mean.s2hashkeepers=apply(s2$hasheskeepers,1,mean)
mean.s3hashkeepers=apply(s3$hasheskeepers,1,mean)
mean.s4hashkeepers=apply(s4$hasheskeepers,1,mean)
mean.hashkeepers.proper=(mean.s1hashkeepers[which(names(mean.s1hashkeepers)=="TS.1")]+
                         mean.s2hashkeepers[which(names(mean.s2hashkeepers)=="TS.2")]+
                         mean.s3hashkeepers[which(names(mean.s3hashkeepers)=="TS.3")]+
                         mean.s4hashkeepers[which(names(mean.s4hashkeepers)=="TS.4")])/4

var.s1hashkeepers=apply(s1$hasheskeepers,1,var)
var.s2hashkeepers=apply(s2$hasheskeepers,1,var)
var.s3hashkeepers=apply(s3$hasheskeepers,1,var)
var.s4hashkeepers=apply(s4$hasheskeepers,1,var)
var.hashkeepers.proper=(var.s1hashkeepers[which(names(var.s1hashkeepers)=="TS.1")]+
                         var.s2hashkeepers[which(names(var.s2hashkeepers)=="TS.2")]+
                         var.s3hashkeepers[which(names(var.s3hashkeepers)=="TS.3")]+
                         var.s4hashkeepers[which(names(var.s4hashkeepers)=="TS.4")])/4

mean.hashkeepers.improper=(mean.s1hashkeepers[which(names(mean.s1hashkeepers)!="TS.1")]+
                         mean.s2hashkeepers[which(names(mean.s2hashkeepers)!="TS.2")]+
                         mean.s3hashkeepers[which(names(mean.s3hashkeepers)!="TS.3")]+
                         mean.s4hashkeepers[which(names(mean.s4hashkeepers)!="TS.4")])/4

var.hashkeepers.improper=(var.s1hashkeepers[which(names(var.s1hashkeepers)!="TS.1")]+
                          var.s2hashkeepers[which(names(var.s2hashkeepers)!="TS.2")]+
                          var.s3hashkeepers[which(names(var.s3hashkeepers)!="TS.3")]+
                          var.s4hashkeepers[which(names(var.s4hashkeepers)!="TS.4")])/4

###Make harmonized data with the same variance

varharmonized.s1hasheskeepers=rbind(sqrt(var.hashkeepers.proper/var.s1hashkeepers[names(var.s1hashkeepers)=="TS.1"])*s1$hasheskeepers[which(names(var.s1hashkeepers)=="TS.1"),],sqrt(var.hashkeepers.improper/var.s1hashkeepers[names(var.s1hashkeepers)!="TS.1"])*s1$hasheskeepers[which(names(var.s1hashkeepers)!="TS.1"),])
varharmonized.s2hasheskeepers=rbind(sqrt(var.hashkeepers.proper/var.s2hashkeepers[names(var.s2hashkeepers)=="TS.2"])*s2$hasheskeepers[which(names(var.s2hashkeepers)=="TS.2"),],sqrt(var.hashkeepers.improper/var.s2hashkeepers[names(var.s2hashkeepers)!="TS.2"])*s2$hasheskeepers[which(names(var.s2hashkeepers)!="TS.2"),])
varharmonized.s3hasheskeepers=rbind(sqrt(var.hashkeepers.proper/var.s3hashkeepers[names(var.s3hashkeepers)=="TS.3"])*s3$hasheskeepers[which(names(var.s3hashkeepers)=="TS.3"),],sqrt(var.hashkeepers.improper/var.s3hashkeepers[names(var.s3hashkeepers)!="TS.3"])*s3$hasheskeepers[which(names(var.s3hashkeepers)!="TS.3"),])
varharmonized.s4hasheskeepers=rbind(sqrt(var.hashkeepers.improper/var.s4hashkeepers[names(var.s4hashkeepers)!="TS.4"])*s4$hasheskeepers[which(names(var.s4hashkeepers)!="TS.4"),],sqrt(var.hashkeepers.proper/var.s4hashkeepers[names(var.s4hashkeepers)=="TS.4"])*s4$hasheskeepers[which(names(var.s4hashkeepers)=="TS.4"),])

###Make same variance data have the same mean

meanharmonized.varharmonized.s1hasheskeepers=rbind(varharmonized.s1hasheskeepers[names(var.s1hashkeepers)=="TS.1",]+mean.hashkeepers.proper-mean(varharmonized.s1hasheskeepers[names(var.s1hashkeepers)=="TS.1",]),varharmonized.s1hasheskeepers[names(var.s1hashkeepers)!="TS.1",]+mean.hashkeepers.improper-mean(varharmonized.s1hasheskeepers[names(var.s1hashkeepers)!="TS.1",]))
meanharmonized.varharmonized.s2hasheskeepers=rbind(varharmonized.s2hasheskeepers[names(var.s2hashkeepers)=="TS.2",]+mean.hashkeepers.proper-mean(varharmonized.s2hasheskeepers[names(var.s2hashkeepers)=="TS.2",]),varharmonized.s2hasheskeepers[names(var.s2hashkeepers)!="TS.2",]+mean.hashkeepers.improper-mean(varharmonized.s2hasheskeepers[names(var.s2hashkeepers)!="TS.2",]))
meanharmonized.varharmonized.s3hasheskeepers=rbind(varharmonized.s3hasheskeepers[names(var.s3hashkeepers)=="TS.3",]+mean.hashkeepers.proper-mean(varharmonized.s3hasheskeepers[names(var.s3hashkeepers)=="TS.3",]),varharmonized.s3hasheskeepers[names(var.s3hashkeepers)!="TS.3",]+mean.hashkeepers.improper-mean(varharmonized.s3hasheskeepers[names(var.s3hashkeepers)!="TS.3",]))
meanharmonized.varharmonized.s4hasheskeepers=rbind(varharmonized.s4hasheskeepers[names(var.s4hashkeepers)!="TS.4",]+mean.hashkeepers.improper-mean(varharmonized.s4hasheskeepers[names(var.s4hashkeepers)!="TS.4",]),varharmonized.s4hasheskeepers[names(var.s4hashkeepers)=="TS.4",]+mean.hashkeepers.proper-mean(varharmonized.s4hasheskeepers[names(var.s4hashkeepers)=="TS.4",]))

###Now add to objects

s1$hashesharmonized=meanharmonized.varharmonized.s1hasheskeepers
s2$hashesharmonized=meanharmonized.varharmonized.s2hasheskeepers
s3$hashesharmonized=meanharmonized.varharmonized.s3hasheskeepers
s4$hashesharmonized=meanharmonized.varharmonized.s4hasheskeepers

###Now adjust for doublet distribution, mean multipled by 2/3, variance multiplied by 4/3

s1.doublet=s1$hashesharmonized
s2.doublet=s2$hashesharmonized
s3.doublet=s3$hashesharmonized
s4.doublet=s4$hashesharmonized

s1.doublet[names(var.s1hashkeepers)=="TS.1",]=sqrt(4/3)*s1$hashesharmonized[names(var.s1hashkeepers)=="TS.1",]
s2.doublet[names(var.s2hashkeepers)=="TS.2",]=sqrt(4/3)*s2$hashesharmonized[names(var.s2hashkeepers)=="TS.2",]
s3.doublet[names(var.s3hashkeepers)=="TS.3",]=sqrt(4/3)*s3$hashesharmonized[names(var.s3hashkeepers)=="TS.3",]
s4.doublet[names(var.s4hashkeepers)=="TS.4",]=sqrt(4/3)*s4$hashesharmonized[names(var.s4hashkeepers)=="TS.4",]

s1.doublet[names(var.s1hashkeepers)=="TS.1",]=s1.doublet[names(var.s1hashkeepers)=="TS.1",]-mean(s1$hashesharmonized[names(var.s1hashkeepers)=="TS.1",])+(2/3)*mean.hashkeepers.proper
s2.doublet[names(var.s2hashkeepers)=="TS.2",]=s2.doublet[names(var.s2hashkeepers)=="TS.2",]-mean(s2$hashesharmonized[names(var.s2hashkeepers)=="TS.2",])+(2/3)*mean.hashkeepers.proper
s3.doublet[names(var.s3hashkeepers)=="TS.3",]=s3.doublet[names(var.s3hashkeepers)=="TS.3",]-mean(s3$hashesharmonized[names(var.s3hashkeepers)=="TS.3",])+(2/3)*mean.hashkeepers.proper
s4.doublet[names(var.s4hashkeepers)=="TS.4",]=s4.doublet[names(var.s4hashkeepers)=="TS.4",]-mean(s4$hashesharmonized[names(var.s4hashkeepers)=="TS.4",])+(2/3)*mean.hashkeepers.proper

s1$doublet=s1.doublet
s2$doublet=s2.doublet
s3$doublet=s3.doublet
s4$doublet=s4.doublet

###Reverse s4 because second hash is proper

s1$properhashesharmonized=s1$hashesharmonized
s2$properhashesharmonized=s2$hashesharmonized
s3$properhashesharmonized=s3$hashesharmonized
s4$properhashesharmonized=rbind(s4$hashesharmonized[2,],s4$hashesharmonized[1,])

rownames(s1$properhashesharmonized)=c("TS.1","TS.2")
rownames(s2$properhashesharmonized)=c("TS.2","TS.3")
rownames(s3$properhashesharmonized)=c("TS.3","TS.4")
rownames(s4$properhashesharmonized)=c("TS.4","TS.1")

s1$properdoublet=s1$doublet
s2$properdoublet=s2$doublet
s3$properdoublet=s3$doublet
s4$properdoublet=rbind(s4$doublet[2,],s4$doublet[1,])

rownames(s1$properdoublet)=c("TS.1","TS.2")
rownames(s2$properdoublet)=c("TS.2","TS.3")
rownames(s3$properdoublet)=c("TS.3","TS.4")
rownames(s4$properdoublet)=c("TS.4","TS.1")

save(s1,s2,s3,s4,file="../data/data2SimulateFrom.RData")
