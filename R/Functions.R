# Modelo
# Alpha=exp(M)
# Eta=exp(U)
# M~N(X'PsiM,vM)
# vM~GI(am,bm)
# Theta~N(m0,c0)
# beta~a*Beta(nu*tau,(1-tau)*nu), a>0
# tau~Beta(ar,br)

sintonizar=function(taxa,tau,mat,i){

  mat=as.matrix(mat)



  mater=(1/50)*sum(mat)

  if(mater>=taxa){
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)-delta
    temp5=exp(temp4)
    return(temp5)
  }else{
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)+delta
    temp5=exp(temp4)
    return(temp5)
  }


}
#
sintonizarN=function(taxa,tau,mat,i){

  mat=as.matrix(mat)



  mater=(1/50)*sum(mat)

  if(mater>=taxa){
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)+delta
    temp5=exp(temp4)
    return(temp5)
  }else{
    delta=min(0.01,(i/50+1)^(-0.5))
    temp4=log(tau)-delta
    temp5=exp(temp4)
    return(temp5)
  }


}
#
#
# b escalar
# def matriz nx2

gCorr<-function(b,def){
  n<-nrow(def)
  R<-exp(-b*(as.matrix(dist(def))))
  mat<-R
  mat


}

# b escalar
# v escalar
# def matriz nx2

gSigma<-function(b,v,def){
  n<-nrow(def)
  R<-exp(-b*(as.matrix(dist(def))))
  mat<-v*R
  mat


}
#
# vec vector
# beta>0
#
R=function(vec,beta){

  if(length(vec)==1){
    return(0)
  }else{

    temp=sum(exp(-beta*(vec[length(vec)]-vec[1:(length(vec)-1)]) ))
    temp
  }

}
#
# alpha>0
# beta>0
# beta>alpha
# lambda>0
# eta>0
logveroWeibull=function(alpha,beta,lambda,eta,tempos){

  sum1=0
  for(i in 1:length(tempos)){
    temp1=log(lambda*eta*tempos[i]^(eta-1)+alpha*R(tempos[1:i],beta))
    sum1=temp1+sum1
  }
  sum2=length(tempos)-R(tempos,beta)
  res=sum1-lambda*tempos[length(tempos)]^(eta)-(alpha/beta)*sum(sum2)
  res
}
#
# M matriz nx1
# beta matriz nx1
# W matriz nx1
# U matriz nx1
# yT matriz mxn
# NN matriz nx1
# ddelta>0 escalar
# 0<f<2pi
# theta>0
LoglikehoodHawkesWeibull=function(M,beta,W,U,yT,NN){

  n<-ncol(yT)
  res<-NULL
  for(i in 1:n){
    tem<-logveroWeibull(exp(M[i,]),beta[i,],exp(W[i,]),exp(U[i,]),yT[1:NN[i,1],i])
    res<-c(res,tem)
  }
  res

}
#
# Amostrar W
# W=log(delta) matriz nx1
# Beta  matriz nx1
# M=log(Alpha) matriz nx1
# exp(M)<Beta
# U=log(delta) matriz nx1
# loca matriz nx2
# XX matriz nxp
# PPs matriz px1
# vv escalar
# bb escalar
# yT matrix mxn
# TT matriz nx1
# NN matriz nx1
samplerW_Weibull<-function(W,beta,M,U,loca,XX,PPs,bb,vv,NN,yT,TT,ff){

  n<-nrow(W)
  WWprop<-as.matrix(mvrnorm(1,W,ff*diag(1,n)))

  SSig<-gSigma(bb,vv,loca)

  postWW<-sum(LoglikehoodHawkesWeibull(M,beta,W,U,yT,NN))-0.5*t(W-XX%*%PPs)%*%solve(SSig)%*%(W-XX%*%PPs)
  postWWprop=sum(LoglikehoodHawkesWeibull(M,beta,WWprop,U,yT,NN))-0.5*t(WWprop-XX%*%PPs)%*%solve(SSig)%*%(WWprop-XX%*%PPs)

  prob=min(exp((postWWprop)-(postWW)),1)


  u=runif(1,0,1)

  if(u<prob){

    Wprox=WWprop

    rejei=1


  }else{

    Wprox=W
    rejei=0
  }





  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res

}
#
# Amostrar U
# W=log(delta) matriz nx1
# Beta  matriz nx1
# M=log(Alpha) matriz nx1
# exp(M)<Beta
# U=exp(eta) matriz nx1
# loca matriz nx2
# XX matriz nxp
# PPs matriz px1
# vv escalar
# bb escalar
# yT matrix mxn
# TT matriz nx1
# NN matriz nx1

samplerU_Weibull<-function(W,beta,M,U,loca,XX,PPs,bb,vv,NN,yT,ff,TT){

  n<-nrow(U)
  Uprop<-as.matrix(mvrnorm(1,U,ff*diag(1,n)))

  SSig<-gSigma(bb,vv,loca)

  postWW<-sum(LoglikehoodHawkesWeibull(M,beta,W,U,yT,NN))-0.5*t(U-XX%*%PPs)%*%solve(SSig)%*%(U-XX%*%PPs)
  postWWprop=sum(LoglikehoodHawkesWeibull(M,beta,W,Uprop,yT,NN))-0.5*t(Uprop-XX%*%PPs)%*%solve(SSig)%*%(Uprop-XX%*%PPs)

  prob=min(exp((postWWprop)-(postWW)),1)


  u=runif(1,0,1)

  if(u<prob){

    Wprox=Uprop

    rejei=1


  }else{

    Wprox=U
    rejei=0
  }





  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res

}

#
# Amostrar M
# W=log(delta) matriz nx1
# Beta  matriz nx1
# M=log(Alpha) matriz nx1
# U=log(Eta)
# exp(M)<Beta
# loca matriz nx2
# XX matriz nxp
# PPs matriz px1
# vv escalar
# bb escalar
# yT matrix mxn
# TT matriz nx1
# NN matriz nx1

samplerM_Weibull<-function(W,beta,M,U,loca,XX,PPs,bb,vv,NN,yT,ff,TT){

  n<-nrow(M)
  Mprop<-as.matrix(mvrnorm(1,M,ff*diag(1,n)))

  if(sum( ifelse(Mprop>exp(beta),1,0) )>=1){
    res=list(M,0)
    return(res)
  }else{

  }

  SSig<-gSigma(bb,vv,loca)

  postWW<-sum(LoglikehoodHawkesWeibull(M,beta,W,U,yT,NN))-0.5*t(M-XX%*%PPs)%*%solve(SSig)%*%(M-XX%*%PPs)
  postWWprop=sum(LoglikehoodHawkesWeibull(Mprop,beta,W,U,yT,NN))-0.5*t(Mprop-XX%*%PPs)%*%solve(SSig)%*%(Mprop-XX%*%PPs)

  prob=min(exp((postWWprop)-(postWW)),1)


  u=runif(1,0,1)

  if(u<prob){

    Wprox=Mprop

    rejei=1


  }else{

    Wprox=M
    rejei=0
  }





  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res

}
#
#
#
samplerbwu=function(WU,v,b,loca,ab,bb,X,Psi,u1){

  bprop=rgamma(1,shape=b*u1, rate = u1)

  SSigprop=gSigma(bprop,v,loca)

  if((det(SSigprop)==0)|(bprop< 0.005)){
    return(list(b,0))
  }

  SSig=gSigma(b,v,loca)
  SSigprop=gSigma(bprop,v,loca)


  logp=-0.5*t(WU-X%*%Psi)%*%solve(SSig)%*%(WU-X%*%Psi)-0.5*log(det(SSig))+(ab-1)*log(b)-bb*b

  logpprop=-0.5*t(WU-X%*%Psi)%*%solve(SSigprop)%*%(WU-X%*%Psi)-0.5*log(det(SSigprop))+(ab-1)*log(bprop)-bb*bprop

  logprob=logpprop+log(dgamma(b,shape=bprop*u1,rate=u1))-(logp+log(dgamma(bprop,shape=b*u1,rate=u1)))
  prob<-min(c(1,exp(logprob)))

  u=runif(1,0,1)

  if(u<prob){

    bprox=bprop

    rejei=1


  }else{

    bprox=b

    rejei=0;

  }



  res=list(bprox,rejei)
  res
}
# Amostrar Beta
# Lambda matriz escalar
# Beta  matriz escalar
# M=log(Alpha) escalar
# M<Q
# XX vector nx1
# PPs matriz px1
# yT matrix vector
Samplerbeta_Weibull<-function(Lambda,Beta,M,U,tau,nu,yT,ff,A,nn){

  Betaprop<-A*rbeta(1, (Beta/A)*ff, (1-(Beta/A))*ff)

  if((M>exp(Betaprop))||(round(Betaprop,6))==A){
    res=list(Beta,0)
    return(res)
  }else{

  }

  postMM<-logveroWeibull(exp(M),Beta,Lambda,exp(U),yT[1:nn])+dbeta((Beta/A),tau*nu,(1-tau)*nu,log = T)
  postMMprop<-logveroWeibull(exp(M),Betaprop,Lambda,exp(U),yT[1:nn])+dbeta((Betaprop/A),tau*nu,(1-tau)*nu,log = T)

  logprob=postMMprop+(1/A)*dbeta( (Beta/A), (Betaprop/A)*ff, (1-(Betaprop/A))*ff,log=T )-(postMM +(1/A)*dbeta( (Betaprop/A), (Beta/A)*ff, (1-(Beta/A))*ff,log=T ) )

  prob<-min(c(1,exp(logprob)))


  u=runif(1,0,1)

  if(u<prob){

    Wprox=Betaprop

    rejei=1


  }else{

    Wprox=Beta
    rejei=0
  }


  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res

}
#
SampleTau=function(tau,nu,Beta,A,atau,btau,ff){

  tauprop<-rbeta(1,tau*ff,(1-tau)*ff)

  if((round(tauprop,6))==1){
    res=list(tau,0)
    return(res)
  }else{

  }

  postMM<-sum(dbeta((Beta/A),tau*nu,(1-tau)*nu,log = T))+dbeta(tau,atau,btau,log=T)
  postMMprop<-sum(dbeta((Beta/A),tauprop*nu,(1-tauprop)*nu,log = T))+dbeta(tauprop,atau,btau,log=T)

  logprob=postMMprop+dbeta(tau, tauprop*ff, (1-tauprop)*ff,log=T )-(postMM +dbeta( tauprop, tau*ff, (1-tau)*ff,log=T ) )

  prob<-min(c(1,exp(logprob)))


  u=runif(1,0,1)

  if(u<prob){

    Wprox=tauprop

    rejei=1


  }else{

    Wprox=tau
    rejei=0
  }


  res=as.matrix(Wprox)
  res=list(Wprox,rejei)
  res


}


################# interpolate mean ########
#
Tempos=function(tt,tat){

  if(tt<tat[1]){

    return(tt)

  }

  cont=1
  output=NULL


  while((tt>tat[cont])&(cont<length(tat))){

    output=c(output,tat[cont])
    cont=cont+1
  }

  output=c(output,tt)
  output
}
#

FunW=function(MatD,p){

  output=(1/MatD^p)/sum(1/MatD^p)
  output

}

#

MeanFunction=function(gam,eta,alpha,beta,t,tant){


  output=gam*t^eta+(alpha/beta)*((length(tant)-R(tant,beta)))
  output


}
#

