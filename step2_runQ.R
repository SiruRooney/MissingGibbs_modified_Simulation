#step2
setwd("F://Missingdatacodes//MissingGibbs2")
rm(list=ls())

load("F://Missingdatacodes//MissingGibbs2//step1_misgibsimdata2.RData")

source("F://Missingdatacodes//MissingGibbs2//Gibbssampler.R")

n.complete=NROW(na.omit(t(cbind(zz,xx))))#18 covariates without missing values
#ind.zx=matrix(0,2+2*NROW(NCOL(cbind(zz,xx))-n.complete),2*NCOL(cbind(zz,xx,y)))
#x1=apply(x,1,function(x){list(x[x!=0])})


ind.zx=list(rep(1,NCOL(cbind(zz,xx))+1),rep(1,n.complete+1),rep(1,n.complete+2),rep(1,n.complete+4),rep(1,n.complete+5),rep(1,n.complete+6))


ptm <- proc.time();ptm
tab.out=Gibbs.sampler(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),cbind(rr.x,ss.y),ind.zx,35)
save.image("F://Missingdatacodes//MissingGibbs2//step2_runQ.RData")
proc.time()-ptm