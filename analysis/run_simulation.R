###Code to actually run simulation

genotype.conversion=function(x,y)
{
    nx=length(x)
    if(nx!=length(y)) stop ("x and y different lengths")
    output=rep(NA,nx)
    output[which(x==0 & y==0)]=0
    output[which(x==0 & y==1)]=1
    output[which(x==1 & y==0)]=1
    output[which(x==1 & y==1)]=1
    output[which(x==1 & is.na(y))]=1
    output[which(is.na(x) & y==1)]=1
    output[which(x==0 & is.na(y))]=0
    output[which(is.na(x) & y==0)]=0
    output[which(is.na(x) & is.na(y))]=NA
    output
}

sim=function(datasets.list,n.cells,prop.doublets,improper.multiplier)
{
    n.datasets=length(datasets.list)
###    print(n.datasets)
    n.data.cells=rep(NA,n.datasets)
    for(j in 1:n.datasets) n.data.cells[j]=ncol(datasets.list[[j]]$hasheskeepers)
###    print(n.data.cells)
    cumsum.cells=cumsum(n.data.cells)
    sum.cells=cumsum.cells[n.datasets]
###    print(cumsum.cells)
    n.doublets=round(prop.doublets*n.cells)
###Sample cells
    samp=sample(1:sum.cells,n.cells,replace=TRUE)
###Initialize cells to sample 1
    which.sample=rep(1,n.cells)
    for(j in 2:n.datasets)
    {
###If cells corresponds to sample j!=1, update the sample id
        which.j=which(samp>cumsum.cells[j-1] & samp<=cumsum.cells[j])
        which.sample[which.j]=j
###Update the cell to be positioned correctly
        samp[which.j]=samp[which.j]-cumsum.cells[j-1]
    }
###    print(which.sample)
###Now sample doublets
    doublet.samp=rep(NA,n.cells)
    doublet.which.sample=rep(NA,n.cells)
###Loop over doublets
    if(n.doublets>0)
    {
        for(k in 1:n.doublets)
        {
###            print(k)
###            print(which.sample[k])
            keep.sampling=TRUE
            while(keep.sampling)
        {
            new.cell=sample(1:sum.cells,1)
###        print(new.cell)
###                print(new.cell)
            which.j=1
            for(j in 2:n.datasets)
            {
                if(new.cell>cumsum.cells[j-1] & new.cell<=cumsum.cells[j])
                {
                    which.j=j
                    break
                }
            }
###        print(which.j)
            if(which.j!=which.sample[k])
            {
                doublet.which.sample[k]=which.j
                doublet.samp[k]=new.cell
                if(which.j>1) doublet.samp[k]=new.cell-cumsum.cells[j-1]
                keep.sampling=FALSE
            }
        }
        }
    }

###    print(doublet.which.sample)
###Do genotypes
    genotype.output=matrix(NA,nrow(datasets.list[[j]]$mutkeepers),n.cells)
    depthtotal.output=matrix(NA,nrow(datasets.list[[j]]$mutkeepers),n.cells)
    depthalt.output=matrix(NA,nrow(datasets.list[[j]]$mutkeepers),n.cells)
    for(k in 1:n.cells)
    {
###            print(k)
###            print(which.sample[k])
###            print(length(datasets.list[[which.sample[k]]]$mutkeepershashkeepers[,samp[k]]))
        genotype.output[,k]=datasets.list[[which.sample[k]]]$mutkeepershashkeepers[,samp[k]]
        depthtotal.output[,k]=datasets.list[[which.sample[k]]]$depthtotalvariantkeepershashkeepers[,samp[k]]
        depthalt.output[,k]=datasets.list[[which.sample[k]]]$depthaltvariantkeepershashkeepers[,samp[k]]
    }
    final.genotype.output=genotype.output
    final.depthtotal.output=depthtotal.output
    final.depthalt.output=depthalt.output
    doublet.output=NULL
    doublet.depthtotal.output=NULL
    doublet.depthalt.output=NULL
###Do doublet genotypes
###Update final genotypes using conversion function
    if(n.doublets>0)
    {
        doublet.output=matrix(NA,nrow(datasets.list[[j]]$mutkeepers),n.doublets)
        doublet.depthtotal.output=matrix(NA,nrow(datasets.list[[j]]$mutkeepers),n.doublets)
        doublet.depthalt.output=matrix(NA,nrow(datasets.list[[j]]$mutkeepers),n.doublets)
###        print(dim(doublet.output))
        for(k in 1:n.doublets)
        {
###            print("k")
###            print(k)
###            print(doublet.samp[k])
###            print(doublet.which.sample[k])
###            print(length(doublet.which.sample))
###            print(datasets.list[[doublet.which.sample[k]]]$mutkeepershashkeepers[,doublet.samp[k]])
            doublet.output[,k]=datasets.list[[doublet.which.sample[k]]]$mutkeepershashkeepers[,doublet.samp[k]]
            doublet.depthtotal.output[,k]=datasets.list[[doublet.which.sample[k]]]$depthtotalvariantkeepershashkeepers[,doublet.samp[k]]
            doublet.depthalt.output[,k]=datasets.list[[doublet.which.sample[k]]]$depthaltvariantkeepershashkeepers[,doublet.samp[k]]
            final.genotype.output[,k]=genotype.conversion(genotype.output[,k],doublet.output[,k])
            final.depthtotal.output[,k]=(depthtotal.output[,k]+doublet.depthtotal.output[,k])/2
            final.depthalt.output[,k]=(depthalt.output[,k]+doublet.depthalt.output[,k])/2
        }
        }
###Do hashes
    hash.output=matrix(NA,n.datasets,n.cells)
    if(n.doublets>0)
        for(k in 1:n.doublets)
        {
            new.sample1=which.sample[k]
###            print(new.sample1)
            new.sample2=doublet.which.sample[k]
###            print(new.sample2)
            new.cell1=samp[k]
###            print(new.cell1)
            new.cell2=doublet.samp[k]
###            print(new.cell2)
###            print(datasets.list[[new.sample1]]$properdoublets[1,new.cell1])
###Fill with improper hashes, then write over with proper hashes in two spots
            hash.output[,k]=datasets.list[[new.sample1]]$properhashesharmonized[2,new.cell1]
            hash.output[new.sample1,k]=datasets.list[[new.sample1]]$properdoublet[1,new.cell1]
            hash.output[new.sample2,k]=datasets.list[[new.sample2]]$properdoublet[1,new.cell2]
###Add a little noise to the improper
            if(n.datasets==3)
            {
                which.improper=(1:3)[-sort(c(new.sample1,new.sample2))]
                hash.output[which.improper,k]=hash.output[which.improper,k]+rnorm(n=1,mean=0,sd=sqrt(improper.multiplier*abs(hash.output[which.improper,k])))
            }
            else if(n.datasets==4)
            {
                which.improper=(1:4)[-sort(c(new.sample1,new.sample2))]
                for(j in 1:2) hash.output[which.improper[j],k]=hash.output[which.improper[j],k]+rnorm(n=1,mean=0,sd=sqrt(improper.multiplier*abs(hash.output[which.improper[j],k])))
            }

        }
        for(k in (n.doublets+1):n.cells)
        {
            new.sample1=which.sample[k]
            new.cell1=samp[k]
###Fill with improper hashes, then write over with proper hashes in one spot
            hash.output[,k]=datasets.list[[new.sample1]]$properhashesharmonized[2,new.cell1]
            hash.output[new.sample1,k]=datasets.list[[new.sample1]]$properhashesharmonized[1,new.cell1]
            if(n.datasets==2)
            {
                which.improper=(1:2)[-new.sample1]
                hash.output[which.improper,k]=hash.output[which.improper,k]+rnorm(n=1,mean=0,sd=sqrt(improper.multiplier*abs(hash.output[which.improper,k])))
            }
            else if(n.datasets==3)
            {
                which.improper=(1:3)[-new.sample1]
                for(j in 1:2) hash.output[which.improper[j],k]=hash.output[which.improper[j],k]+rnorm(n=1,mean=0,sd=sqrt(improper.multiplier*abs(hash.output[which.improper[j],k])))
            }
            else if(n.datasets==4)
            {
                which.improper=(1:4)[-new.sample1]
                for(j in 1:3) hash.output[which.improper[j],k]=hash.output[which.improper[j],k]+rnorm(n=1,mean=0,sd=sqrt(improper.multiplier*abs(hash.output[which.improper[j],k])))
            }
        }
    #truth=which.sample
    truth=as.character(which.sample)
    if(n.doublets > 0)
    {
        for(k in 1:n.doublets)
        {
            samples=sort(c(which.sample[k],doublet.which.sample[k]))
###            print(samples)
            truth[k]= paste0(samples[1],",",samples[2])
        }
    }
    list(samp=samp,which.sample=which.sample,doublet.samp=doublet.samp,doublet.which.sample=doublet.which.sample,n.doublets=n.doublets,final.genotype.output=final.genotype.output,genotype.output=genotype.output,doublet.output=doublet.output,hash.output=hash.output,final.depthtotal.output=final.depthtotal.output,depthtotal.output=depthtotal.output,doublet.depthtotal.output=doublet.depthtotal.output,final.depthalt.output=final.depthalt.output,depthalt.output=depthalt.output,doublet.depthalt.output=doublet.depthalt.output,truth=truth)
}
###samp is the cell number (vector of length n.cells)
###which.sample is the sample number (vector of length n.cells)
###doublet.samp is the doublet cell number (vector of length n.cells with NAs for non-doublets)
###doublet.which.samp is the doublet sample (vector of length n.cells with NAs for non-doublets)
###n.doublets is the number of doublets (a single value)
###final.genotype.output is the genotypes including revision for doublets (matrix variants by n.cells)
###genotype.output is the genotypes just for the initial sample so not including doublets (matrix variants by n.cells)
###doublet.output is the genotypes for the doublet (matrix variants by n.doublets)
###final.depthtotal.output is the depthtotals including revision for doublets (matrix variants by n.cells)
###depthtotal.output is the depthtotals just for the initial sample so not including doublets (matrix variants by n.cells)
###doublet.depthtotal.output is the depthtotals for the doublet (matrix variants by n.doublets)
###final.depthtotal.output is the depthtotals including revision for doublets (matrix variants by n.cells)
###depthtotal.output is the depthtotals just for the initial sample so not including doublets (matrix variants by n.cells)
###doublet.depthalt.output is the depthtotals for the doublet (matrix variants by n.doublets)
###hash.output is the hash output (matrix n.datasets by n.cells)
###truth (vector of length n.cells)

####################################################################
####################################################################
