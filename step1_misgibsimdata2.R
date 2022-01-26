setwd("F://Missingdatacodes//MissingGibbs2")
rm(list=ls())


zz=rnorm(600*15)
zz=matrix(zz,nrow=600,ncol=15)
xx=rbinom(600*5,1,.5)
xx=matrix(xx,nrow=600,ncol=5)

Alpha1=c(1,-1,1.8);Alpha2=c(1,1,-0.9,-1.5)
Beta=c(1,1,1,-1,2,-1.7);Phi1=c(1,1,.78,-2,1.8);Phi2=c(1,-1,1.5,1.5,-1.85)
Eta=c(1,1,1,-1,1.6,-2)

for (i in 1:NROW(zz)){
  xx[i,3]=rbinom(1,1,exp(c(1,zz[i,c(3,6)])%*%Alpha1)/(1+exp(c(1,zz[i,c(3,6)])%*%Alpha1)))
  xx[i,5]=rbinom(1,1,exp(c(1,zz[i,c(9,10)],xx[i,3])%*%Alpha2)/(1+exp(c(1,zz[i,c(9,10)],xx[i,3])%*%Alpha2)))
}
summary(glm(xx[,3]~zz,family=binomial(link=logit)))
summary(glm(xx[,5]~zz+xx[,-5],family=binomial(link=logit)))

yy=matrix(0,NROW(zz),1)
for (i in 1:NROW(zz)){
  yy[i]=rbinom(1,1,exp(c(1,zz[i,c(1,3,6)],xx[i,c(3,5)])%*%Beta)/(1+exp(c(1,zz[i,c(1,3,6)],xx[i,c(3,5)])%*%Beta)))
}
summary(glm(yy~zz+xx,family=binomial(link=logit)))

rr.x=matrix(0,NROW(zz),2)
for (i in 1:NROW(zz)){
  rr.x[i,1]=rbinom(1,1,exp(c(1,zz[i,c(1,2)],xx[i,3],yy[i])%*%Phi1)/(1+exp(c(1,zz[i,c(1,2)],xx[i,3],yy[i])%*%Phi1)))
  rr.x[i,2]=rbinom(1,1,exp(c(1,zz[i,9],xx[i,c(2,3)],yy[i])%*%Phi2)/(1+exp(c(1,zz[i,9],xx[i,c(2,3)],yy[i])%*%Phi2)))
}
summary(glm(rr.x[,1]~zz+xx+yy,family=binomial(link=logit)))
summary(glm(rr.x[,2]~zz+xx+yy+rr.x[,1],family=binomial(link=logit)))

ss.y=matrix(0,NROW(zz),1)
for(i in 1:NROW(zz)){
  ss.y[i]=rbinom(1,1,exp(c(1,zz[i,c(1:3)],xx[i,2],yy[i])%*%Eta)/(1+exp(c(1,zz[i,1:3],xx[i,2],yy[i])%*%Eta)))
}
summary(glm(ss.y~zz+xx+yy+rr.x,family=binomial(link=logit)))

xx[,3][(rr.x[,1]==0)]=NA;xx[,5][(rr.x[,2]==0)]=NA
yy[(ss.y==0)]=NA

save.image("F://Missingdatacodes//MissingGibbs2//step1_misgibsimdata2.RData")