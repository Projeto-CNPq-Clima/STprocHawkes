#' Estimate Parameters of Spatio-Temporal Hawkes Model using MCMC
#'
#' @description
#' This function estimates the parameters of a spatio-temporal Hawkes model
#' using Markov Chain Monte Carlo (MCMC) methods.
#'
#' @param data A vector of event times for the phenomenon of interest. This input is similar to the `data` argument in the `STModelWeibullMCMC` function.
#' @param sites A matrix of geographic coordinates corresponding to the event locations. This input is similar to the `sites` argument in the `STModelWeibullMCMC` function.
#' @param Xw A matrix of covariates associated with the \(W\) process. Defaults to a matrix combining an intercept and `sites`.
#' @param Xm A matrix of covariates associated with the \(M\) process. Defaults to the same as `Xw`.
#' @param Xu A matrix of covariates associated with the \(U\) process. Defaults to the same as `Xw`.
#' @param prior A list of hyperparameters defining the prior distributions in the model.
#' @param iteration Number of MCMC iterations. Defaults to 1000.
#' @param burnin Number of iterations for the burn-in period. Defaults to 900.
#' @param jump Thinning parameter for the MCMC sampling. Defaults to 10.
#'
#' @return A list containing:
#'   - Estimated parameters for the spatio-temporal Hawkes model.
#'   - MCMC samples for each parameter.
#'   - Convergence diagnostics and summary statistics.
#'
#' @export
STModelHawkesMCMC<-function(data, sites,Xw=as.matrix(cbind(as.matrix(rep(1,n)),sites)),Xm=Xw,Xu=Xw,
                            prior=list(Nu=2,Aw=as.matrix(rep(0,ncol(Xw))),
                                       Bw=diag(100,ncol(Xw)),
                                       Am=as.matrix(rep(0,ncol(Xm))),
                                       Bm=diag(100,ncol(Xm)),
                                       Au=as.matrix(rep(0,ncol(Xu))),
                                       Bu=diag(100,ncol(Xu)),
                                       Aa=2,
                                       aaw=0.1,
                                       bbw=0.1,
                                       aam=1,
                                       bbm=1,
                                       aau=1,
                                       bbu=1,
                                       aa1w=0.1,
                                       bb1w=0.1,
                                       aa1m=0.1,
                                       bb1m=0.1,
                                       aa1u=0.1,
                                       bb1u=0.1),iteration,burnin,jump){
  n=ncol(data)
  m=nrow(data)
  #  Xw=as.matrix(cbind(as.matrix(rep(1,n)),sites))
  #  Xm=Xw
  #  Xu=Xw
  tempdados=is.na(data)
  nj=as.matrix(m-apply(tempdados,2,sum))



  Tt=array(NA,dim=c(1,n))

  for(y in 1:n){
    Tt[1,y]=data[nj[y],y]
  }

  Tt=t(Tt)

  # Parametros iniciais

  Ww=as.matrix(seq(1,1.5, length=n),ncol(data))
  Mm=as.matrix(seq(0.01,0.02, length=n),ncol(data))
  Beta=as.matrix(seq(1.1,1.2, length=n),ncol(data))
  Uu=as.matrix(rep(-1,ncol(data)))
  Tau=0.5
  Nu=2
  Psiw=as.matrix(rep(0,ncol(Xw)))
  Psim=as.matrix(rep(0,ncol(Xm)))
  Psiu=as.matrix(rep(0,ncol(Xu)))

  bw=0.1
  bm=0.1
  bu=0.1
  vw=0.1
  vm=0.1
  vu=0.1
  # Hiper
  Nu=prior$Nu
  Aw=prior$Aw
  Bw=prior$Bw
  Am=prior$Am
  Bm=prior$Bm
  Au=prior$Au
  Bu=prior$Bu

  Aa=prior$Aa

  aaw=prior$aaw
  bbw=prior$bbw
  aam=prior$aam
  bbm=prior$bbm
  aau=prior$aau
  bbu=prior$bbu


  aa1w=prior$aa1w
  bb1w=prior$bb1w
  aa1m=prior$aa1m
  bb1m=prior$bb1m
  aa1u=prior$aa1u
  bb1u=prior$bb1u

  # iteration=420000
  #  burnin=400000
  #  jump=10

  SU1=1/1000
  SU2=1/1000
  SU3=1000
  SU4=1000
  SU5=rep(2000,n)
  SU6=1/1000
  SU7=1000
  MW=NULL
  MWT=NULL
  MU=NULL
  MUT=NULL

  MPsiw=NULL
  MMT=NULL
  MbwT=NULL
  Mbw=NULL
  MbmT=NULL
  Mbm=NULL
  Mbu=NULL
  MbuT=NULL

  MPsim=NULL
  Mvw=NULL
  Mvm=NULL
  Mvu=NULL

  MBeta=NULL
  MTau=NULL
  MBetaT=NULL
  MM=NULL
  MPsiu=NULL


  for(i in 1:iteration){

    if(i<=burnin){

      temp=tryCatch(samplerW_Weibull(Ww,Beta,Mm,Uu,sites,Xw,Psiw,bw,vw,nj,data,t(Tt),SU1),error = function(e) e)
      if(is.matrix(temp[[1]])==T){
        Ww=temp[[1]]
        MWT=c(MWT,temp[[2]])
      }else{
        MWT=c(MWT,0)
      }

      varPsiw=solve(solve(Bw)+t(Xw)%*%solve(gSigma(bw,vw,sites))%*%Xw)
      medPsiw=(t(Aw)%*%solve(Bw)+t(Ww)%*%solve(gSigma(bw,vw,sites))%*%Xw)%*%varPsiw
      Psiw=as.matrix(MASS::mvrnorm(1,medPsiw,varPsiw))

      RRw=gCorr(bw,sites)
      aaaw=(n/2)+aaw
      bbbw=0.5*t(Ww-Xw%*%Psiw)%*%solve(RRw)%*%(Ww-Xw%*%Psiw)+bbw
      vw=1/rgamma(1,shape=aaaw, rate = bbbw)

      temp=samplerbwu(Ww,vw,bw,sites,aa1w,bb1w,Xw,Psiw,SU3)
      bw=temp[[1]]
      MbwT=c(MbwT,temp[[2]])


      temp=tryCatch(samplerM_Weibull(Ww,Beta,Mm,Uu,sites,Xm,Psim,bm,vm,nj,data,SU2,t(Tt)),error = function(e) e)
      if(is.matrix(temp[[1]])==T){
        Mm=temp[[1]]
        MMT=c(MMT,temp[[2]])
      }else{
        MMT=c(MMT,0)
      }

      varPsim=solve(solve(Bm)+t(Xm)%*%solve(gSigma(bm,vm,sites))%*%Xm)
      medPsim=(t(Am)%*%solve(Bm)+t(Mm)%*%solve(gSigma(bm,vm,sites))%*%Xm)%*%varPsim
      Psim=as.matrix(MASS::mvrnorm(1,medPsim,varPsim))

      RRm=gCorr(bm,sites)
      aaam=(n/2)+aam
      bbbm=0.5*t(Mm-Xm%*%Psim)%*%solve(RRm)%*%(Mm-Xm%*%Psim)+bbm
      vm=1/rgamma(1,shape=aaam, rate = bbbm)

      temp=samplerbwu(Mm,vm,bm,sites,aa1m,bb1m,Xm,Psim,SU4)
      bm=temp[[1]]
      MbmT=c(MbmT,temp[[2]])


      temp=tryCatch(samplerU_Weibull(Ww,Beta,Mm,Uu,sites,Xu,Psiu,bu,vu,nj,data,SU6,t(Tt)),error = function(e) e)
      if(is.matrix(temp[[1]])==T){
        Uu=temp[[1]]
        MUT=c(MUT,temp[[2]])
      }else{
        MUT=c(MUT,0)
      }

      varPsiu=solve(solve(Bu)+t(Xu)%*%solve(gSigma(bu,vu,sites))%*%Xu)
      medPsiu=(t(Au)%*%solve(Bu)+t(Uu)%*%solve(gSigma(bu,vu,sites))%*%Xu)%*%varPsiu
      Psiu=as.matrix(MASS::mvrnorm(1,medPsiu,varPsiu))

      RRu=gCorr(bu,sites)
      aaau=(n/2)+aau
      bbbu=0.5*t(Uu-Xu%*%Psiu)%*%solve(RRu)%*%(Uu-Xu%*%Psiu)+bbu
      vu=1/rgamma(1,shape=aaau, rate = bbbu)

      temp=samplerbwu(Uu,vu,bu,sites,aa1u,bb1u,Xu,Psiu,SU7)
      bu=temp[[1]]
      MbuT=c(MbuT,temp[[2]])


      auxbetaT=NULL
      for(j in 1:n){
        temp=tryCatch(Samplerbeta_Weibull(exp(Ww[j,1]),Beta[j,1],Mm[j,1],Uu[j,1],Tau,Nu,data[,j],SU5[j],Aa,nj[j,1]),error = function(e) e)
        if(is.numeric(temp[[1]])==T){
          Beta[j,1]=temp[[1]]
          auxbetaT=c(auxbetaT,temp[[2]])
        }else{
          auxbetaT=c(auxbetaT,0)
        }
      }

      MBetaT=rbind(MBetaT,t(as.matrix(auxbetaT)))

      temp=SampleTau(Tau,Nu,Beta,Aa,0.5,0.5,Aa)
      Tau=temp[[1]]

      if((i%%50)==0){
        SU1=sintonizarN(0.25,SU1,MWT[(i-50+1):i],i)
        SU2=sintonizarN(0.25,SU2,MMT[(i-50+1):i],i)
        SU3=sintonizar(0.4,SU3,MbwT[(i-50+1):i],i)
        SU4=sintonizar(0.4,SU4,MbmT[(i-50+1):i],i)
        SU6=sintonizarN(0.25,SU6,MUT[(i-50+1):i],i)
        SU7=sintonizar(0.4,SU7,MbuT[(i-50+1):i],i)

      }else{

      }




    }else{


      if((i%%jump)==0){

        temp=tryCatch(samplerW_Weibull(Ww,Beta,Mm,Uu,sites,Xw,Psiw,bw,vw,nj,data,t(Tt),SU1),error = function(e) e)
        if(is.matrix(temp[[1]])==T){
          Ww=temp[[1]]
          MW=rbind(MW,t(Ww))
        }else{
          MW=rbind(MW,t(Ww))
        }

        varPsiw=solve(solve(Bw)+t(Xw)%*%solve(gSigma(bw,vw,sites))%*%Xw)
        medPsiw=(t(Aw)%*%solve(Bw)+t(Ww)%*%solve(gSigma(bw,vw,sites))%*%Xw)%*%varPsiw
        Psiw=as.matrix(MASS::mvrnorm(1,medPsiw,varPsiw))
        MPsiw=rbind(MPsiw,t(Psiw))

        RRw=gCorr(bw,sites)
        aaaw=(n/2)+aaw
        bbbw=0.5*t(Ww-Xw%*%Psiw)%*%solve(RRw)%*%(Ww-Xw%*%Psiw)+bbw
        vw=1/rgamma(1,shape=aaaw, rate = bbbw)
        Mvw=c(Mvw,vw)

        temp=samplerbwu(Ww,vw,bw,sites,aa1w,bb1w,Xw,Psiw,SU3)
        bw=temp[[1]]
        Mbw=c(Mbw,bw)
        temp=tryCatch(samplerM_Weibull(Ww,Beta,Mm,Uu,sites,Xm,Psim,bm,vm,nj,data,SU2,t(Tt)),error = function(e) e)
        if(is.matrix(temp[[1]])==T){
          Mm=temp[[1]]
          MM=rbind(MM,t(Mm))
        }else{
          MM=rbind(MM,t(Mm))
        }

        varPsim=solve(solve(Bm)+t(Xm)%*%solve(gSigma(bm,vm,sites))%*%Xm)
        medPsim=(t(Am)%*%solve(Bm)+t(Mm)%*%solve(gSigma(bm,vm,sites))%*%Xm)%*%varPsim
        Psim=as.matrix(MASS::mvrnorm(1,medPsim,varPsim))
        MPsim=rbind(MPsim,t(Psim))

        RRm=gCorr(bm,sites)
        aaam=(n/2)+aam
        bbbm=0.5*t(Mm-Xm%*%Psim)%*%solve(RRm)%*%(Mm-Xm%*%Psim)+bbm
        vm=1/rgamma(1,shape=aaam, rate = bbbm)
        Mvm=c(Mvm,vm)

        temp=samplerbwu(Mm,vm,bm,sites,aa1m,bb1m,Xm,Psim,SU4)
        bm=temp[[1]]
        Mbm=c(Mbm,bm)

        temp=tryCatch(samplerU_Weibull(Ww,Beta,Mm,Uu,sites,Xu,Psiu,bu,vu,nj,data,SU6,t(Tt)),error = function(e) e)
        if(is.matrix(temp[[1]])==T){
          Uu=temp[[1]]
          MU=rbind(MU,t(Uu))
        }else{
          MU=rbind(MU,t(Uu))
        }

        varPsiu=solve(solve(Bu)+t(Xu)%*%solve(gSigma(bu,vu,sites))%*%Xu)
        medPsiu=(t(Au)%*%solve(Bu)+t(Uu)%*%solve(gSigma(bu,vu,sites))%*%Xu)%*%varPsiu
        Psiu=as.matrix(MASS::mvrnorm(1,medPsiu,varPsiu))
        MPsiu=rbind(MPsiu,t(Psiu))

        RRu=gCorr(bu,sites)
        aaau=(n/2)+aau
        bbbu=0.5*t(Uu-Xu%*%Psiu)%*%solve(RRu)%*%(Uu-Xu%*%Psiu)+bbu
        vu=1/rgamma(1,shape=aaau, rate = bbbu)
        Mvu=c(Mvu,vu)

        temp=samplerbwu(Uu,vu,bu,sites,aa1u,bb1u,Xu,Psiu,SU7)
        bu=temp[[1]]
        Mbu=c(Mbu,bu)

        for(j in 1:n){
          temp=tryCatch(Samplerbeta_Weibull(exp(Ww[j,1]),Beta[j,1],Mm[j,1],Uu[j,1],Tau,Nu,data[,j],SU5[j],Aa,nj[j,1]),error = function(e) e)
          if(is.numeric(temp[[1]])==T){
            Beta[j,1]=temp[[1]]
          }else{
          }
        }
        MBeta=rbind(MBeta,t(Beta))

        temp=SampleTau(Tau,Nu,Beta,Aa,0.5,0.5,Aa)
        Tau=temp[[1]]
        MTau=c(MTau,Tau)

      }else{
        temp=tryCatch(samplerW_Weibull(Ww,Beta,Mm,Uu,sites,Xw,Psiw,bw,vw,nj,data,t(Tt),SU1),error = function(e) e)
        if(is.matrix(temp[[1]])==T){
          Ww=temp[[1]]
        }else{
        }

        varPsiw=solve(solve(Bw)+t(Xw)%*%solve(gSigma(bw,vw,sites))%*%Xw)
        medPsiw=(t(Aw)%*%solve(Bw)+t(Ww)%*%solve(gSigma(bw,vw,sites))%*%Xw)%*%varPsiw
        Psiw=as.matrix(MASS::mvrnorm(1,medPsiw,varPsiw))

        RRw=gCorr(bw,sites)
        aaaw=(n/2)+aaw
        bbbw=0.5*t(Ww-Xw%*%Psiw)%*%solve(RRw)%*%(Ww-Xw%*%Psiw)+bbw
        vw=1/rgamma(1,shape=aaaw, rate = bbbw)

        temp=samplerbwu(Ww,vw,bw,sites,aa1w,bb1w,Xw,Psiw,SU3)
        bw=temp[[1]]

        temp=tryCatch(samplerM_Weibull(Ww,Beta,Mm,Uu,sites,Xm,Psim,bm,vm,nj,data,SU2,t(Tt)),error = function(e) e)
        if(is.matrix(temp[[1]])==T){
          Mm=temp[[1]]
        }else{
        }

        varPsim=solve(solve(Bm)+t(Xm)%*%solve(gSigma(bm,vm,sites))%*%Xm)
        medPsim=(t(Am)%*%solve(Bm)+t(Mm)%*%solve(gSigma(bm,vm,sites))%*%Xm)%*%varPsim
        Psim=as.matrix(MASS::mvrnorm(1,medPsim,varPsim))

        RRm=gCorr(bm,sites)
        aaam=(n/2)+aam
        bbbm=0.5*t(Mm-Xm%*%Psim)%*%solve(RRm)%*%(Mm-Xm%*%Psim)+bbm
        vm=1/rgamma(1,shape=aaam, rate = bbbm)

        temp=samplerbwu(Mm,vm,bm,sites,aa1m,bb1m,Xm,Psim,SU4)
        bm=temp[[1]]

        temp=tryCatch(samplerU_Weibull(Ww,Beta,Mm,Uu,sites,Xu,Psiu,bu,vu,nj,data,SU6,t(Tt)),error = function(e) e)
        if(is.matrix(temp[[1]])==T){
          Uu=temp[[1]]
        }else{
        }

        varPsiu=solve(solve(Bu)+t(Xu)%*%solve(gSigma(bu,vu,sites))%*%Xu)
        medPsiu=(t(Au)%*%solve(Bu)+t(Uu)%*%solve(gSigma(bu,vu,sites))%*%Xu)%*%varPsiu
        Psiu=as.matrix(MASS::mvrnorm(1,medPsiu,varPsiu))

        RRu=gCorr(bu,sites)
        aaau=(n/2)+aau
        bbbu=0.5*t(Uu-Xu%*%Psiu)%*%solve(RRu)%*%(Uu-Xu%*%Psiu)+bbu
        vu=1/rgamma(1,shape=aaau, rate = bbbu)

        temp=samplerbwu(Uu,vu,bu,sites,aa1u,bb1u,Xu,Psiu,SU7)
        bu=temp[[1]]

        for(j in 1:n){
          temp=tryCatch(Samplerbeta_Weibull(exp(Ww[j,1]),Beta[j,1],Mm[j,1],Uu[j,1],Tau,Nu,data[,j],SU5[j],Aa,nj[j,1]),error =function(e) e)
          if(is.numeric(temp[[1]])==T){
            Beta[j,1]=temp[[1]]
          }else{
          }
        }

        temp=SampleTau(Tau,Nu,Beta,Aa,0.5,0.5,Aa)
        Tau=temp[[1]]


      }



    }


    print(i)

  }

  resul<-list(MW, MWT, MU, MUT,MPsiw,MMT,MbwT,Mbw, MbmT,Mbm,Mbu,MbuT,MPsim,Mvw,Mvm,Mvu,MBeta,MTau,MBetaT,MM,MPsiu)
  names(resul)<-c("MW", "MWT", "MU", "MUT","MPsiw","MMT","MbwT","Mbw", "MbmT","Mbm","Mbu","MbuT","MPsim","Mvw","Mvm","Mvu","MBeta","MTau","MBetaT","MM","MPsiu")
  return(resul)
}
