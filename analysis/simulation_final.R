###Simulate the minimum number of cells by sample

dat=read.delim("summary_cell_experiment_1-7_14-16.txt",sep="\t",header=TRUE,row.names=1)

###Cell counts: Mean 5300, SD=2500, minimum=1500

pmult=as.numeric(dat[29,]/dat[1,])
nsamps=c(1,1,1,1,2,3,4,4,8,8)
fit=lm(pmult~nsamps)
fit.sigma=summary(fit)$sigma

###Proportion multiplets: .1+0.025*nsamps

p2=c(733/(1577+733),1577/(1577+733))
p3=c(570/(570+767+1049),767/(570+767+1049),1049/(570+767+1049))
p4=c(567/(4313+567+2598+1402),1402/(4313+567+2598+1402),2598/(4313+567+2598+1402),4313/(4313+567+2598+1402),
     1140/(1140+3325+2748+2414),2414/(1140+3325+2748+2414),2748/(1140+3325+2748+2414),3325/(1140+3325+2748+2414))
p8=c(228/(524+213+378+228+446+231+315+2768+2399),
     231/(524+213+378+228+446+231+315+2768+2399),
     315/(524+213+378+228+446+231+315+2768+2399),
     378/(524+213+378+228+446+231+315+2768+2399),
     446/(524+213+378+228+446+231+315+2768+2399),
     524/(524+213+378+228+446+231+315+2768+2399),
     2399/(524+213+378+228+446+231+315+2768+2399),
     2768/(524+213+378+228+446+231+315+2768+2399),
     45/(463+275+354+45+234+435+142+605),
     142/(463+275+354+45+234+435+142+605),
     234/(463+275+354+45+234+435+142+605),
     275/(463+275+354+45+234+435+142+605),
     354/(463+275+354+45+234+435+142+605),
     435/(463+275+354+45+234+435+142+605),
     463/(463+275+354+45+234+435+142+605),
     605/(463+275+354+45+234+435+142+605))

pmins=c(NA,NA,NA,NA,p2[1],p3[1],p4[1],p4[5],p8[1],p8[9])
fit2=lm(pmins~nsamps)
fit2.sigma=summary(fit2)$sigma

###Proportion of cells in smallest count sample: .33-0.04nsamps

simfunc=function(n.sims=10000,n.samps,mean1,sd1,min1,sd2,sd3,min3)
{
    n.cells=round(rnorm(n.sims,mean=mean1,sd=sd1))
    n.cells[which(n.cells<min1)]=min1
    percentage.multiplet=rnorm(n.sims,mean=.1+0.025*n.samps,sd=sd2)
    n.multiplet=round(percentage.multiplet*n.cells)
    final.cells=n.cells-n.multiplet
    final.cells
    percentage.minimum=rnorm(n.sims,mean=.33-0.04*n.samps,sd=sd3)
    percentage.minimum[which(percentage.minimum<min3)]=min3
    minimum.cells=round(percentage.minimum*final.cells)
    minimum.cells
}

set.seed(12345)
simfunc.2=simfunc(n.sims=10000,n.samps=2,mean1=5300,sd1=2500,min1=1500,sd2=.067,sd3=0.070,min3=0.02)
simfunc.3=simfunc(n.sims=10000,n.samps=3,mean1=5300,sd1=2500,min1=1500,sd2=.067,sd3=0.070,min3=0.02)
simfunc.4=simfunc(n.sims=10000,n.samps=4,mean1=5300,sd1=2500,min1=1500,sd2=.067,sd3=0.070,min3=0.02)
simfunc.5=simfunc(n.sims=10000,n.samps=5,mean1=5300,sd1=2500,min1=1500,sd2=.067,sd3=0.070,min3=0.02)
simfunc.6=simfunc(n.sims=10000,n.samps=6,mean1=5300,sd1=2500,min1=1500,sd2=.067,sd3=0.070,min3=0.02)
simfunc.7=simfunc(n.sims=10000,n.samps=7,mean1=5300,sd1=2500,min1=1500,sd2=.067,sd3=0.070,min3=0.02)
simfunc.8=simfunc(n.sims=10000,n.samps=8,mean1=5300,sd1=2500,min1=1500,sd2=.067,sd3=0.070,min3=0.02)

p.100=c(sum(simfunc.2>=100)/10000,sum(simfunc.3>=100)/10000,sum(simfunc.4>=100)/10000,sum(simfunc.5>=100)/10000,sum(simfunc.6>=100)/10000,sum(simfunc.7>=100)/10000,sum(simfunc.8>=100)/10000)
p.200=c(sum(simfunc.2>=200)/10000,sum(simfunc.3>=200)/10000,sum(simfunc.4>=200)/10000,sum(simfunc.5>=200)/10000,sum(simfunc.6>=200)/10000,sum(simfunc.7>=200)/10000,sum(simfunc.8>=200)/10000)
p.300=c(sum(simfunc.2>=300)/10000,sum(simfunc.3>=300)/10000,sum(simfunc.4>=300)/10000,sum(simfunc.5>=300)/10000,sum(simfunc.6>=300)/10000,sum(simfunc.7>=300)/10000,sum(simfunc.8>=300)/10000)
p.400=c(sum(simfunc.2>=400)/10000,sum(simfunc.3>=400)/10000,sum(simfunc.4>=400)/10000,sum(simfunc.5>=400)/10000,sum(simfunc.6>=400)/10000,sum(simfunc.7>=400)/10000,sum(simfunc.8>=400)/10000)
p.500=c(sum(simfunc.2>=500)/10000,sum(simfunc.3>=500)/10000,sum(simfunc.4>=500)/10000,sum(simfunc.5>=500)/10000,sum(simfunc.6>=500)/10000,sum(simfunc.7>=500)/10000,sum(simfunc.8>=500)/10000)

pdf("simplot_final.pdf")

plot(c(2,8),c(0,1),type="n",xlab="Number of Samples Multiplexed",ylab="Estimated Probability of Minimum Number of Cells")
points(2:8,p.100,pch=16,col=1)
lines(2:8,p.100,col=1)
points(2:8,p.200,pch=16,col=2)
lines(2:8,p.200,col=2)
points(2:8,p.300,pch=16,col=3)
lines(2:8,p.300,col=3)
points(2:8,p.400,pch=16,col=4)
lines(2:8,p.400,col=4)
points(2:8,p.500,pch=16,col=5)
lines(2:8,p.500,col=5)

legend(2.5,0.4,legend=c(">=100 Cells",">=200 Cells",">=300 Cells",">=400 Cells",">=500 Cells"),lty=rep(1,5),col=1:5)

dev.off()
