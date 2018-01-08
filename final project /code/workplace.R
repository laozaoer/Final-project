###########
generateData=function(Mu,Sigma,beta,sig,n1,n2,n3){
  p=length(beta)
  Xtrain=mvrnorm(n1,Mu,Sigma)
  ytrain=Xtrain%*%beta+rnorm(n1,0,sig)
  Xvalid=mvrnorm(n2,Mu,Sigma)
  yvalid=Xvalid%*%beta+rnorm(n2,0,sig)
  Xtest=mvrnorm(n3,rep(0,p),Sigma)
  ytest=Xtest%*%beta+rnorm(n3,0,sig)
  return(list(Xtrain=Xtrain,ytrain=ytrain,Xvalid=Xvalid,yvalid=yvalid,Xtest=Xtest,ytest=ytest))
}

#####################
simulate=function(Mu,Sigma,beta,sig,n1,n2,n3,nsim,Lambda,a){
  p=length(beta)
  softBetahat=matrix(rep(0,nsim*p),ncol=p)
  softError=rep(0,nsim)
  hardBetahat=matrix(rep(0,nsim*p),ncol=p)
  hardError=rep(0,nsim)
  scadBetahat=matrix(rep(0,nsim*p),ncol=p)
  scadError=rep(0,nsim)
  
  for(i in 1:nsim){
    Result=generateData(Mu,Sigma,beta,sig,n1,n2,n3)
    Xtrain=Result$Xtrain
    ytrain=Result$ytrain
    Xvalid=Result$Xvalid
    yvalid=Result$yvalid
    Xtest=Result$Xtest
    ytest=Result$ytest
    
    lambda=softTISPcrossValidation(Xvalid,yvalid,5,Lambda)
    softBetahat[i,]=softTISP(Xtrain,ytrain,lambda)
    yhat=Xtest%*%softBetahat[i,]
    error=sum((ytest-yhat)*(ytest-yhat))/n3 #MSE
    softError[i]=error
    
    lambda=hardTISPcrossValidation(Xvalid,yvalid,5,Lambda)
    hardBetahat[i,]=hardTISP(Xtrain,ytrain,lambda)
    yhat=Xtest%*%hardBetahat[i,]
    error=sum((ytest-yhat)*(ytest-yhat))/n3 #MSE
    hardError[i]=error
    
    lambda=scadTISPcrossValidation(Xvalid,yvalid,5,Lambda,a)
    scadBetahat[i,]=scadTISP(Xtrain,ytrain,lambda,a)
    yhat=Xtest%*%scadBetahat[i,]
    error=sum((ytest-yhat)*(ytest-yhat))/n3 #MSE
    scadError[i]=error
  }
  return(list(softBetahat=softBetahat,softMSE=softError,
              hardBetahat=hardBetahat,hardMSE=hardError,
              scadBetahat=scadBetahat,scadMSE=scadError))
}


#####################
sparseError=function(Betahat,beta){
  p=length(beta)
  nsim=nrow(Betahat)
  k=0
  for(j in 1:p){
    for(i in 1:nsim){
      if(sign(Betahat[i,j])==sign(beta[j]))
        k=k+1
    }
  }
  sparseError=1-k/(nsim*p)
  return(sparseError)
}
#########################
properZero=function(Betahat,beta){
  p=length(beta)
  nsim=nrow(Betahat)
  zeroLoca=rep(0,p)
  for(i in 1:p){
    if(beta[i]==0)
      zeroLoca[i]=i
  }
  zeroLoca=zeroLoca[zeroLoca>0]
  zeroNum=length(zeroLoca)
  k=0
  for(j in 1:zeroNum){
    for(i in 1:nsim){
      if(Betahat[i,zeroLoca[j]]==0)
        k=k+1
    }
  }
  properZero=k/(nsim*zeroNum)
  return(properZero)
}
#############################
properNonzero=function(Betahat,beta){
  p=length(beta)
  nsim=nrow(Betahat)
  nonzeroLoca=rep(0,p)
  for(i in 1:p){
    if(beta[i]!=0)
      nonzeroLoca[i]=i
  }
  nonzeroLoca=nonzeroLoca[nonzeroLoca>0]
  nonzeroNum=length(nonzeroLoca)
  k=0
  for(j in 1:nonzeroNum){
    for(i in 1:nsim){
      if(Betahat[i,nonzeroLoca[j]]!=0)
        k=k+1
    }
  }
  properNonzero=k/(nsim*nonzeroNum)
  return(properNonzero)
}


##########################parameters setting
beta=c(3,1.5,0,0,2,0,0,0) #c(3,1.5,0,0,2,rep(0,95)) p=100
p=length(beta)
sig=2 #8
n1=20;n2=100;n3=200
rho=0.5 #0.85      ##### 4 kinds of settings{(sig=2,rho=0.5),...(sig=8,rho=0.85)}

covariance=function(rho,p){
  Sigma=matrix(rep(0,p^2),ncol=p)
  for(i in 1:p){
    for(j in i:p){
      Sigma[i,j]=rho^(abs(i-j))
      Sigma[j,i]=Sigma[i,j]
    }
  }
  return(Sigma)
}
Sigma=covariance(rho,p)
Mu=rep(0,p)
library(MASS)
Lambda=seq(0,5,0.2)
nsim=1000
a=3.7 #scad

###############################Simulation
Outcome=simulate(Mu,Sigma,beta,sig,n1,n2,n3,nsim,Lambda,a)

softBetahat=Outcome$softBetahat
softMSE=Outcome$softMSE
hardBetahat=Outcome$hardBetahat
hardMSE=Outcome$hardMSE
scadBetahat=Outcome$scadBetahat
scadMSE=Outcome$scadMSE

##############performance
performComp=matrix(rep(0,12),ncol=3,dimnames = list(c("Mse", "Sparse error","proper zeros",
                                                      "proper nonzeros"),
                                                    c("softTISP(Lasso)", "hardTISP", 
                                                      "scadTISP")))
    performComp[1,1]=mean(softMSE)
    performComp[2,1]=sparseError(softBetahat,beta)
    performComp[3,1]=properZero(softBetahat,beta)
    performComp[4,1]=properNonzero(softBetahat,beta)
    
    performComp[1,2]=mean(softMSE)
    performComp[2,2]=sparseError(hardBetahat,beta)
    performComp[3,2]=properZero(hardBetahat,beta)
    performComp[4,2]=properNonzero(hardBetahat,beta)
    
    performComp[1,3]=mean(scadMSE)
    performComp[2,3]=sparseError(scadBetahat,beta)
    performComp[3,3]=properZero(scadBetahat,beta)
    performComp[4,3]=properNonzero(scadBetahat,beta)

    
################# make tables: when p<n (p=8,n=20) nsim=1000
write.csv(softBetahat,file="文档/statistical compution/final/Lasso estimations(sig=2,rho=0.5)")
write.csv(performComp,file="文档/statistical compution/final/Performance comparisons(p=8,sig=2,
          rho=0.5)")
write.csv(performComp,file="文档/statistical compution/final/Performance comparisons(p=8,sig=8,
          rho=0.5)")
write.csv(performComp,file="文档/statistical compution/final/Performance comparisons(p=8,sig=2,
          rho=0.85)")
write.csv(performComp,file="文档/statistical compution/final/Performance comparisons(p=8,sig=8,
          rho=0.85)")


##################### make tables: when p>n (p=100,n=20) nsim=100
write.csv(performComp,file="文档/statistical compution/final/Performance comparisons(p=100,
sig=2, rho=0.5)")
write.csv(performComp,file="文档/statistical compution/final/Performance comparisons(p=100,
sig=8,rho=0.5)")
write.csv(performComp,file="文档/statistical compution/final/Performance comparisons(p=100,
sig=2,rho=0.85)")
write.csv(performComp,file="文档/statistical compution/final/Performance comparisons(p=100,
sig=8,rho=0.85)")

