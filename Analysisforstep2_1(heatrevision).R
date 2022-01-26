#load("F:/Missingdatacodes/MissingGibbs2_modified/MissingGibbs2/step2_runQ2.RData")
#tab.obs=tab.out
#rm(tab.out)
load("F:/Missingdatacodes/MissingGibbs2_modified/MissingGibbs2/step2_runQ3.RData")
tab.com=tab.out

for (i.tab in 1:6){tab.com[[1]][[i.tab]]=rbind(tab.com[[1]][[i.tab]],tab.out3[[1]][[i.tab]])}
rm(tab.out,tab.out3)

#load("F:/Missingdatacodes/MissingGibbs2_modified/MissingGibbs2/step2_runcomp.RData")

load("F:/Missingdatacodes/MissingGibbs2_modified/MissingGibbs2/step2_runcomp2.RData")
tab.obs=tab.out
for (i.tab in 1:6){tab.obs[[1]][[i.tab]]=rbind(tab.obs[[1]][[i.tab]],tab.out2[[1]][[i.tab]])}
rm(tab.out,tab.out2)

nrow.tab=7380
#par(mfrow=c(2,1))

plot(tab.com[[1]][[1]][1:nrow.tab,1],type="l",col="blue",ylim=c(3150,3850),ylab="BIC_{Q}",xlab="iteration")#
abline(h=min(tab.com[[1]][[1]][1:nrow.tab,1]),col="red")

axis(2,seq(3100,3900,100),seq(3100,3900,100))
title(main="I-chart for the generated BIC sequences")
par(new=T)
plot(tab.obs[[1]][[1]][1:nrow.tab,1],type="l",col="darkgreen",ylim=c(850,1000),ylab="",xlab="",xaxt="n",yaxt="n")
abline(h=min(tab.obs[[1]][[1]][1:nrow.tab,1]),col="orange")
axis(4,seq(850,1000,50),seq(850,1000,50))
mtext("BIC_{obs}",4)
legend("topright",lty=c(1,1),col=c("blue","darkgreen"),legend=c("BIC_{Q}","BIC_{obs}"))

###############################################################################
table(tab.com[[1]][[1]][1:nrow.tab,1])
NROW(table(tab.com[[1]][[1]][1:nrow.tab,1]))#511
table(tab.obs[[1]][[1]][1:nrow.tab,1])
NROW(table(tab.obs[[1]][[1]][1:nrow.tab,1]))#470
####################################################################

BAPE=function(list.bape){
  Beta=rep(0,sum(list.bape[[1]]))
  Eta=rep(0,sum(list.bape[[NROW(list.bape)]]))
  Alpha=lapply(list.bape[2:(1+0.5*(NROW(list.bape)-2))],function(x){y=rep(0,sum(x));return(y)})
  Phi=lapply(list.bape[(2+0.5*(NROW(list.bape)-2)):(NROW(list.bape)-1)],function(x){y=rep(0,sum(x));return(y)})
  return(list(Beta,Alpha,Phi,Eta))
}
ind.BAPE=function(list.bape){
  ind.beta=list.bape[[1]]
  ind.eta=list.bape[[NROW(list.bape)]]
  ind.alpha=lapply(list.bape[2:(1+0.5*(NROW(list.bape)-2))],identity)
  ind.phi=lapply(list.bape[(2+0.5*(NROW(list.bape)-2)):(NROW(list.bape)-1)],identity)
  return(list(ind.beta,ind.alpha,ind.phi,ind.eta))
}

EM.BIC=function(predata,bape,ind.bape,rs.ind){
  #first and second derivative functions for Italydata
  #combination for missing data
  Bitmatrix<-function(n){
    set<-0:(2^n-1)
    rst<-matrix(0,ncol=n,nrow=2^n)
    for (i in 1:n){
      rst[,i]=ifelse((set-rowSums(rst*rep(c(2^((n-1):0)),each=2^n)))/(2^(n-i))>=1,1,0)
    }
    rst
  }
  
  
  ###################################################################################
  #conditional weighted probability
  con.w=function(predamis.k,bape,ind.bape,rs.k){#(z.intpt,x,y,bape,ind.bape,r,s,miscate)
    
    p.ymis=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]])))
    p.ymis=as.vector(p.ymis)
    denom.k=dbinom(predamis.k[,NCOL(predamis.k)],1,p.ymis)
    
    predamisr.k=matrix(0.0,nrow=NROW(predamis.k),ncol=NCOL(predamis.k)+2,byrow=TRUE);
    predamisr.k[,1:NCOL(predamis.k)]=predamis.k;predamisr.k[,-(1:NCOL(predamis.k))]=matrix(rs.k[-NROW(rs.k)],nrow=NROW(predamis.k),ncol=NROW(rs.k)-1,byrow=TRUE)
    for (i.xr in 1:(NROW(rs.k)-1)){
      p.xrmis=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)])%*%as.matrix(bape[[2]][[i.xr]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)])%*%as.matrix(bape[[2]][[i.xr]])))
      p.xrmis=as.vector(p.xrmis)
      denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[2]][[i.xr]])+1],1,p.xrmis)
      
      #print(dim(as.matrix(predamisr.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])))
      #print(i.xr);
      #if(i.xr==2){print(as.matrix(bape[[3]][i.xr]))}else{;}
      p.xrmis=exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])%*%as.matrix(bape[[3]][[i.xr]]))/(1+exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])%*%as.matrix(bape[[3]][[i.xr]])))
      p.xrmis=as.vector(p.xrmis)
      denom.k=denom.k*dbinom(rs.k[i.xr],1,p.xrmis)
    }
    
    p.smis=exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predamisr.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
    p.smis=as.vector(p.smis)
    denom.k=denom.k*dbinom(rs.k[NROW(rs.k)],1,p.smis)
    mistype.p=denom.k/sum(denom.k)
    
    as.vector(mistype.p)
  }
  
  ###############################################################################
  # create the function to generate covatiates list
  
  # the indicator of which sample
  comp.ind=complete.cases(predata)
  mis.type=sapply(c(1:3),Bitmatrix)
  
  
  Dy=function(predata,bape,ind.bape,rs.ind,comp.ind,mis.type){
    predata.inpt=matrix(0.0,NROW(predata),NCOL(predata)+1);
    predata.inpt[,1]=rep(1,NROW(predata));predata.inpt[,-1]=predata
    #observed data
    p.y=exp(as.matrix(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]]))/(1+exp(as.matrix(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]])))
    p.y=as.vector(p.y)
    
    D1y=colSums((predata.inpt[comp.ind==TRUE,NROW(ind.bape[[1]])+1]-p.y)*predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])
    D2y.k=t(apply(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,ind.bape[[1]]==1],1,function(x){x%o%x}))*(-p.y*(1-p.y))
    D2y.k=colSums(D2y.k)
    D2y=matrix(D2y.k,sum(ind.bape[[1]]),sum(ind.bape[[1]]))
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predata.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predata.inpt),byrow=TRUE)
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k,bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
      
      p.y=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)])%*%as.matrix(bape[[1]])))
      p.y=as.vector(p.y)
      D1y=D1y+colSums(as.vector(predamis.k[,NCOL(predamis.k)]-p.y)*predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]*con.P)
      D2y.k=t(apply(predamis.k[,1:NROW(ind.bape[[1]])][,ind.bape[[1]]==1],1,function(x){x%o%x}))*(-p.y)*(1-p.y)*con.P
      D2y.k=colSums(D2y.k)
      D2y=D2y+matrix(D2y.k,sum(ind.bape[[1]]),sum(ind.bape[[1]]))
    }
    list(D1y,D2y)
  }
  
  
  Dx=function(predata,bape,ind.bape,rs.ind,comp.ind,mistype,i.xr){
    predata.inpt=matrix(0.0,NROW(predata),NCOL(predata)+1);
    predata.inpt[,1]=rep(1,NROW(predata));predata.inpt[,-1]=predata  
    D1x=list(rep(0,NROW(bape[[2]][[1]])),rep(0,NROW(bape[[2]][[2]])));
    D2x=list(matrix(0,NROW(bape[[2]][[1]]),NROW(bape[[2]][[1]])),matrix(0,NROW(bape[[2]][[2]]),NROW(bape[[2]][[2]])));
    
    #observed data
    p.x=exp(as.matrix(predata.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])%*%as.matrix(bape[[2]][[i.xr]]))/
      (1+exp(as.matrix(predata.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])%*%as.matrix(bape[[2]][[i.xr]])))
    p.x=as.vector(p.x) 
    D1x[[i.xr]]=colSums((predata.inpt[comp.ind==TRUE,NROW(ind.bape[[2]][[i.xr]])+1]-p.x)*predata.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])
    D2x.k=t(apply(predata.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.x*(1-p.x))
    D2x.k=colSums(D2x.k)
    D2x[[i.xr]]=matrix(D2x.k,sum(ind.bape[[2]][[i.xr]]),sum(ind.bape[[2]][[i.xr]])) 
    
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predata.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predata.inpt),byrow=TRUE)
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k,bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
      
      p.x=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)])%*%as.matrix(bape[[2]][[i.xr]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)])%*%as.matrix(bape[[2]][[i.xr]])))
      p.x=as.vector(p.x)
      D1x[[i.xr]]=D1x[[i.xr]]+colSums(as.vector(predamis.k[,NROW(ind.bape[[2]][[i.xr]])+1]-p.x)*predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]*con.P)
      D2x.k=t(apply(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.x)*(1-p.x)*con.P
      D2x.k=colSums(D2x.k)
      D2x[[i.xr]]=D2x[[i.xr]]+matrix(D2x.k,sum(ind.bape[[2]][[i.xr]]),sum(ind.bape[[2]][[i.xr]]))
    }
    
    list(D1x,D2x)
  }
  
  Dr=function(predata,bape,ind.bape,rs.ind,comp.ind,mistype,i.xr){
    predatar.inpt=matrix(0.0,NROW(predata),NCOL(predata)+NCOL(rs.ind))
    predatar.inpt[,1]=rep(1,NROW(predata));
    predatar.inpt[,2:(NCOL(predata)+1)]=predata;predatar.inpt[,(NCOL(predatar.inpt)-1):(NCOL(predatar.inpt))]=rs.ind[,-NCOL(rs.ind)]
    D1r=list(rep(0,NROW(bape[[3]][[1]])),rep(0,NROW(bape[[3]][[2]])));
    D2r=list(matrix(0,NROW(bape[[3]][[1]]),NROW(bape[[3]][[1]])),matrix(0,NROW(bape[[3]][[2]]),NROW(bape[[3]][[2]])));
    
    #observed data
    p.r=exp(as.matrix(predatar.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])%*%as.matrix(bape[[3]][[i.xr]]))/
      (1+exp(as.matrix(predatar.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])%*%as.matrix(bape[[3]][[i.xr]])))
    p.r=as.vector(p.r) 
    D1r[[i.xr]]=colSums((predatar.inpt[comp.ind==TRUE,NROW(ind.bape[[3]][[i.xr]])+1]-p.r)*predatar.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])
    D2r.k=t(apply(predatar.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.r*(1-p.r))
    D2r.k=colSums(D2r.k)
    D2r[[i.xr]]=matrix(D2r.k,sum(ind.bape[[3]][[i.xr]]),sum(ind.bape[[3]][[i.xr]])) 
    
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){#cat("i.mis",i.mis,"\n");#print(bape[[3]][[2]])
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predatar.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatar.inpt),byrow=TRUE)
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k[,1:(NCOL(predata)+1)],bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
      
      p.r=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])%*%as.matrix(bape[[3]][[i.xr]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)])%*%as.matrix(bape[[3]][[i.xr]])))
      p.r=as.vector(p.r)
      D1r[[i.xr]]=D1r[[i.xr]]+colSums(as.vector(predamis.k[,NROW(ind.bape[[3]][[i.xr]])+1]-p.r)*predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]*con.P)
      D2r.k=t(apply(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1],1,function(x){x%o%x}))*(-p.r)*(1-p.r)*con.P
      D2r.k=colSums(D2r.k)
      D2r[[i.xr]]=D2r[[i.xr]]+matrix(D2r.k,sum(ind.bape[[3]][[i.xr]]),sum(ind.bape[[3]][[i.xr]]))
    }
    
    list(D1r,D2r)
    
  }
  
  Ds=function(predata,bape,ind.bape,rs.ind,comp.ind,mistype){
    predatars.inpt=matrix(0.0,NROW(predata),NCOL(predata)+NCOL(rs.ind)+1)
    predatars.inpt[,1]=rep(1,NROW(predata));
    predatars.inpt[,2:(NCOL(predata)+1)]=predata;predatars.inpt[,(NCOL(predatars.inpt)-2):(NCOL(predatars.inpt))]=rs.ind
    ###################################################################
    #  if(sum(ind.bape[[4]])==1){
    #   D1s=sum((predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[4]])+1]-p.s)*predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])
    #   D2s.k=(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1]^2)*(-p.s*(1-p.s))
    #  D2s.k=sum(D2s.k)
    #  D2s=matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    
    #unobserved data
    #   for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
    #     mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
    #      predamis.k=matrix(predatars.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatars.inpt),byrow=TRUE)
    #      predamis.k[is.na(predamis.k)==T]=mistype.k
    #      con.P=con.w(predamis.k[,1:(NCOL(predata)+1)],bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
    
    #     p.s=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
    #   p.s=as.vector(p.s)
    #   D1s=D1s+sum(as.vector(predamis.k[,NCOL(predamis.k)]-p.s)*predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]*con.P)
    
    #  D2s.k=(predamis.k[,1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1]^2)*(-p.s)*(1-p.s)*con.P
    
    #   D2s.k=sum(D2s.k)
    #   D2s=D2s+matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    # }
    # }else{
    #observed data
    p.s=exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]])/(1+exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]%*%bape[[4]]))
    p.s=as.vector(p.s)
    D1s=colSums((predatars.inpt[comp.ind==TRUE,NROW(ind.bape[[4]])+1]-p.s)*predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])
    D2s.k=t(apply(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1],1,function(x){x%o%x}))*(-p.s*(1-p.s))
    D2s.k=colSums(D2s.k)
    D2s=matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    #unobserved data
    for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){ #cat("i.misds",i.mis,"\n")
      mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
      predamis.k=matrix(predatars.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatars.inpt),byrow=TRUE)
      predamis.k[is.na(predamis.k)==T]=mistype.k
      con.P=con.w(predamis.k[,1:(NCOL(predata)+1)],bape,ind.bape,rs.ind[comp.ind==FALSE,][i.mis,])
      
      p.s=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
      p.s=as.vector(p.s)
      D1s=D1s+colSums(as.vector(predamis.k[,NCOL(predamis.k)]-p.s)*predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)]*con.P)
      D2s.k=t(apply(predamis.k[,1:NROW(ind.bape[[4]])][,ind.bape[[4]]==1],1,function(x){x%o%x}))*(-p.s)*(1-p.s)*con.P
      D2s.k=colSums(D2s.k)
      D2s=D2s+matrix(D2s.k,sum(ind.bape[[4]]),sum(ind.bape[[4]]))
    }
    list(D1s,D2s)
  }
  
  
  ################################################################################################################
  
  i.EM=1
  bape0=bape
  while((i.EM==1)|(max(abs(unlist(bape)-unlist(bape0)))>10e-5)){
    bape0=bape
    cat("i.EM",i.EM,"\n")
    if(sum(ind.bape[[1]])==1){
      bape[[1]]=glm(predata[,NCOL(predata)]~1,family=binomial(link=logit))$coefficients
    }else{
      D12y=Dy(predata,bape,ind.bape,rs.ind,comp.ind,mis.type)
      bape[[1]]=bape[[1]]-solve(D12y[[2]],D12y[[1]])
    }
    
    for (i.mis in 1:(NCOL(rs.ind)-1)){
      if(sum(ind.bape[[2]][[i.mis]])==1){
        bape[[2]][[i.mis]]=glm(predata[,complete.cases(t(predata))==FALSE][,i.mis]~1,family=binomial(link=logit))$coefficients
      }else{
        D12x=Dx(predata,bape,ind.bape,rs.ind,comp.ind,mis.type,i.mis)
        bape[[2]][[i.mis]]=bape[[2]][[i.mis]]-solve(D12x[[2]][[i.mis]],D12x[[1]][[i.mis]])
      }
      
      #print(bape)
      if(sum(ind.bape[[3]][[i.mis]])==1){
        bape[[3]][[i.mis]]=glm(rs.ind[,i.mis]~1,family=binomial(link=logit))$coefficients
      }else{
        D12r=Dr(predata,bape,ind.bape,rs.ind,comp.ind,mis.type,i.mis)
        bape[[3]][[i.mis]]=bape[[3]][[i.mis]]-solve(D12r[[2]][[i.mis]],D12r[[1]][[i.mis]])}
    }
    
    if(sum(ind.bape[[4]])==1){
      bape[[4]]=glm(rs.ind[,3]~1,family=binomial(link=logit))$coefficients
    }else{D12s=Ds(predata,bape,ind.bape,rs.ind,comp.ind,mis.type)
    bape[[4]]=bape[[4]]-solve(D12s[[2]],D12s[[1]])}
    
    i.EM=i.EM+1
    #print(max(abs(unlist(bape)-unlist(bape0))))
  }
  
  cat("i.EM",i.EM,"\n")
  
  # definition of BIC
  BIC.def=function(predata,bape,ind.bape,rs.ind,comp.ind,mis.type){
    p=NROW(unlist(bape))
    #log information of BIC
    log.BIC=function(predata,bape,ind.bape,rs.ind,comp.ind,mis.type){
      predatars.inpt=matrix(0.0,NROW(predata),NCOL(predata)+NCOL(rs.ind)+1)
      predatars.inpt[,1]=rep(1,NROW(predata));
      predatars.inpt[,2:(NCOL(predata)+1)]=predata;predatars.inpt[,(NCOL(predatars.inpt)-2):(NCOL(predatars.inpt))]=rs.ind
      
      #observed data
      p.y=exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]])/(1+exp(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]]))
      p.y=as.vector(p.y)
      L.obs=dbinom(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[1]])+1],1,p.y)
      for (i.xr in 1:2){
        p.x=exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])%*%as.matrix(bape[[2]][[i.xr]]))/
          (1+exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[2]][[i.xr]])][,ind.bape[[2]][[i.xr]]==1])%*%as.matrix(bape[[2]][[i.xr]])))
        p.x=as.vector(p.x)
        L.obs=L.obs*dbinom(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[2]][[i.xr]])+1],1,p.x)
        p.r=exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])%*%as.matrix(bape[[3]][[i.xr]]))/
          (1+exp(as.matrix(predatars.inpt[comp.ind==TRUE,1:NROW(ind.bape[[3]][[i.xr]])][,ind.bape[[3]][[i.xr]]==1])%*%as.matrix(bape[[3]][[i.xr]])))
        p.r=as.vector(p.r)
        L.obs=L.obs*dbinom(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[3]][[i.xr]])+1],1,p.r)
      }
      p.s=exp(as.matrix(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predatars.inpt[(comp.ind==TRUE),1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
      p.s=as.vector(p.s)
      L.obs=L.obs*dbinom(predatars.inpt[(comp.ind==TRUE),NROW(ind.bape[[4]])+1],1,p.s)
      l.obs=sum(log(L.obs))
      obs=l.obs
      #unobserved data
      for (i.mis in 1:NROW(predata[comp.ind==FALSE,])){
        mistype.k=switch(as.character(NCOL(rs.ind)-sum(rs.ind[(comp.ind==FALSE),][i.mis,])),"1"=mis.type[[1]],"2"=mis.type[[2]],"3"=mis.type[[3]])
        predamis.k=matrix(predatars.inpt[(comp.ind==FALSE),][i.mis,],nrow=NROW(mistype.k),ncol=NCOL(predatars.inpt),byrow=TRUE)
        predamis.k[is.na(predamis.k)==T]=mistype.k
        
        p.ymis=exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]])/(1+exp(predamis.k[,1:NROW(ind.bape[[1]])][,(ind.bape[[1]]==1)]%*%bape[[1]]))
        p.ymis=as.vector(p.ymis)
        denom.k=dbinom(predamis.k[,NROW(ind.bape[[1]])+1],1,p.ymis)
        #
        for (i.xr in 1:(NCOL(rs.ind)-1)){
          p.xrmis=exp(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]%*%bape[[2]][[i.xr]])/(1+exp(predamis.k[,1:NROW(ind.bape[[2]][[i.xr]])][,(ind.bape[[2]][[i.xr]]==1)]%*%bape[[2]][[i.xr]]))
          p.xrmis=as.vector(p.xrmis)
          denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[2]][[i.xr]])+1],1,p.xrmis)
          
          p.xrmis=exp(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]%*%bape[[3]][[i.xr]])/(1+exp(predamis.k[,1:NROW(ind.bape[[3]][[i.xr]])][,(ind.bape[[3]][[i.xr]]==1)]%*%bape[[3]][[i.xr]]))
          p.xrmis=as.vector(p.xrmis)
          denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[3]][[i.xr]])+1],1,p.xrmis)
        }
        
        p.smis=exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]]))/(1+exp(as.matrix(predamis.k[,1:NROW(ind.bape[[4]])][,(ind.bape[[4]]==1)])%*%as.matrix(bape[[4]])))
        p.smis=as.vector(p.smis)
        denom.k=denom.k*dbinom(predamis.k[,NROW(ind.bape[[4]])+1],1,p.smis)
        l.obs=l.obs+sum(log(denom.k)*as.vector(denom.k/sum(denom.k)))
      }
      Q=l.obs
      return(list(obs=obs,Q=Q))
    }
    
    lik=log.BIC(predata,bape,ind.bape,rs.ind,comp.ind,mis.type)
    bic.obs=-2*lik$obs+p*log(NROW(predata[comp.ind==TRUE,]))
    bic.Q=-2*lik$Q+p*log(NROW(predata))
    return(list(bic.Q=bic.Q,bic.obs=bic.obs))
  }
  list(BIC.def(predata,bape,ind.bape,rs.ind,comp.ind,mis.type),par.esti=bape)
}
##############################################################################
ind.zx=list(rep(0,NCOL(cbind(zz,xx))+1),rep(0,n.complete+1),rep(0,n.complete+2),rep(0,n.complete+4),rep(0,n.complete+5),rep(0,n.complete+6))
ind.zx[[1]][c(1,2,4,7,20,21)]=1
ind.zx[[2]][c(1,4,7)]=1
ind.zx[[3]][c(1,10,11,20)]=1
ind.zx[[4]][c(1,2,3,20,22)]=1
ind.zx[[5]][c(1,10,18,20,22)]=1
ind.zx[[6]][c(1:4,18,22)]=1

EM.BIC(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),BAPE(ind.zx),ind.BAPE(ind.zx),cbind(rr.x,ss.y))#bic.Q=3353.589,bic.obs=934.5969
#####################################################################################
bic.obs=tab.obs[[1]][[1]][,1]
colSums(tab.obs[[1]][[1]][1:7380,][bic.obs[1:7380]==min(bic.obs[1:7380]),-1])#bic.obs==min(bic.obs)
colSums(tab.obs[[1]][[2]][bic.obs==min(bic.obs),-1])
colSums(tab.obs[[1]][[3]][bic.obs==min(bic.obs),-1])
colSums(tab.obs[[1]][[4]][bic.obs==min(bic.obs),-1])
colSums(tab.obs[[1]][[5]][bic.obs==min(bic.obs),-1])
colSums(tab.obs[[1]][[6]][bic.obs==min(bic.obs),-1])
obs.list=lapply(tab.obs[[1]],function(x,y){colSums(x[y==min(y)[1],-1])},bic.obs)#y==min(y)
obs.m=matrix(0,6,24)
for (i in 1:6){obs.m[i,1:NROW(obs.list[[i]])]=obs.list[[i]]}

bic.com=tab.com[[1]][[1]][,1]
colSums(tab.com[[1]][[1]][bic.com==min(bic.com),-1])
colSums(tab.com[[1]][[2]][bic.com==min(bic.com),-1])
colSums(tab.com[[1]][[3]][bic.com==min(bic.com),-1])
colSums(tab.com[[1]][[4]][bic.com==min(bic.com),-1])
colSums(tab.com[[1]][[5]][bic.com==min(bic.com),-1])
colSums(tab.com[[1]][[6]][bic.com==min(bic.com),-1])
com.list=lapply(tab.com[[1]],function(x,y){colSums(x[y==min(y),-1])},bic.com)
com.m=matrix(0,6,24)
for (i in 1:6){com.m[i,1:NROW(com.list[[i]])]=com.list[[i]]}


true.m=matrix(0,6,24)
for (i in 1:6){true.m[i,1:NROW(ind.zx[[i]])]=ind.zx[[i]]}

#colnames(com.m)=c(" ","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","x1","x2","x3","x4","x5","?","r1","r2")
colnames(com.m)=c(" ","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","x1","x2","y","r1","r2")
rownames(com.m)=c("Q_logit(y)","Q_logit(x4)","Q_logit(x5)","Q_logit(r1)","Q_logit(r2)","Q_logit(r3)")

colnames(obs.m)=c(" ","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","x1","x2","y","r1","r2")
rownames(obs.m)=c("obs_logit(y)","obs_logit(x4)","obs_logit(x5)","obs_logit(r1)","obs_logit(r2)","obs_logit(r3)")

colnames(true.m)=c(" ","z1","z2","z3","z4","z5","z6","z7","z8","z9","z10","z11","z12","z13","z14","z15","z16","z17","z18","x1","x2","y","r1","r2")
rownames(true.m)=c("true_logit(y)","true_logit(x4)","true_logit(x5)","true_logit(r1)","true_logit(r2)","true_logit(r3)")
library(RColorBrewer)
coul = colorRampPalette(brewer.pal(8, "PiYG"))(25)
par(mfrow=c(1,2))
heatmap(t(com.m[,-1]/16),Colv=NA,Rowv=NA,scale="column",col=coul)
#heatmap(t(obs.m[,-1]/259),Colv=NA,Rowv=NA,scale="column",col=coul)
heatmap(t(true.m[,-1]),Colv=NA,Rowv=NA,scale="column",col=coul)
#heatmap(t(rbind(true.m[,-1],rep(0,23),com.m[,-1]/16,rep(0,23),obs.m[,-1]/259)),Colv=NA,Rowv=NA,scale="column",col=coul)
heatmap(t(rbind(true.m[,-1],rep(0,23),com.m[,-1]/16,rep(0,23))),Colv=NA,Rowv=NA,scale="column",col=coul)








########################################################################################################################
source("F:/Missingdatacodes/MissingGibbs2_modified/MissingGibbs2/MisGibparameter.R")
source("F:/Missingdatacodes/MissingGibbs2_modified/MissingGibbs2/fisherinformation.R")
com1.list=lapply(com.list,function(x){x[-1][(x[-1]/16)<0.8]=rep(0,NROW(x[-1][(x[-1]/16)<0.8]));x[-1][(x[-1]/16)>=0.8]=rep(1,NROW(x[-1][(x[-1]/16)>=0.8]));x[1]=1;x})
obs1.list=lapply(obs.list,function(x){x[-1][(x[-1]/259)<0.9]=rep(0,NROW(x[-1][(x[-1]/259)<0.9]));x[-1][(x[-1]/259)>=0.9]=rep(1,NROW(x[-1][(x[-1]/259)>=0.9]));x[1]=1;x})
#(i.EM<=1000)
com.par=EM.BIC1(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),BAPE(com1.list),ind.BAPE(com1.list),cbind(rr.x,ss.y))
#obs.par=EM.BIC1(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),BAPE(obs1.list),ind.BAPE(obs1.list),cbind(rr.x,ss.y))

#comb=list()
#for (i.comb in 1:NROW(com1.list)){comb[[i.comb]]=obs1.list[[i.comb]]+com1.list[[i.comb]]}
#for(i.comb in 1:NROW(com1.list)){comb[[i.comb]][comb[[i.comb]]!=0]=1}
#comb.par=EM.BIC1(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),BAPE(comb),ind.BAPE(comb),cbind(rr.x,ss.y))

com.fish=fisherinfo(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),com.par[[1]],com.par[[2]],cbind(rr.x,ss.y))
#obs.fish=fisherinfo(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),obs.par[[1]],obs.par[[2]],cbind(rr.x,ss.y))
#comb.fish=fisherinfo(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),comb.par[[1]],comb.par[[2]],cbind(rr.x,ss.y))
2-2*pnorm(abs(unlist(comb.par[[1]]))/sqrt(diag(solve(comb.fish))))
sqrt(diag(solve(com.fish)))
2-2*pnorm(abs(unlist(com.par[[1]]))/sqrt(diag(solve(com.fish))))
#the initial full model
ind.zx=list(rep(1,NCOL(cbind(zz,xx))+1),rep(1,n.complete+1),rep(1,n.complete+2),rep(1,n.complete+4),rep(1,n.complete+5),rep(1,n.complete+6))
initial.par=EM.BIC1(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),BAPE(ind.zx),ind.BAPE(ind.zx),cbind(rr.x,ss.y))
initial.fish=fisherinfo(cbind(zz,xx[,c(1,2,4)],xx[,c(3,5)],yy),initial.par[[1]],initial.par[[2]],cbind(rr.x,ss.y))
2-2*pnorm(abs(unlist(initial.par[[1]]))/sqrt(diag(solve(initial.fish))))
