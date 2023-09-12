#######################################################################################################################
## UFO-B: An Utility-based Model Free Phase I/II Dose Optimization Design for Immunotherapy Trials
##          by Yingjie Qiu and Yong Zang
##
## This file contains two functions for simulation studies of UFO designs to find the optimal dose for immunotherapy. 
## get.oc.UFO_B() is a function used to generate operating characteristics for the proposed designs
## get.df.UFO_B() is a function used for dose-finding algorithm of the actual trial.
## get.dose_selection.UFO_B() is a function used for optimal dose selection by the end of the actual trial.
## At the end of this file, examples are provided to illustrate how to use the proposed designs to conduct a clinical trial.
########################################################################################################################


##########################################################################################################
## Function to generate operating characteristics of the UFO-B design 
## To use this function, library "Iso" and 'truncdist' should be installed in R. 
##
##  Arguments:
## pT0.true :true toxicity rates for each dose level without immune response 
## pT1.true :true toxicity rates for each dose level with immune response 
## pE0.true :true toxicity rates for each dose level without immune response 
## pE1.true :true toxicity rates for each dose level with immune response 
## rho0: correlation between efficacy and toxicity for those without immune response
## rho1: correlation between efficacy and toxicity for those with immune response
## pI.true: true immune response rate for each dose level
## targerT: highest acceptable toxicity rate
## targetE: lowest acceptable efficacy rate
## ncohort: total number of cohort
## cohortsize: sample size for each cohort
## startdose: starting dose for dose finding trial
## p.tox: dose deescalation boundary
## cutoff.eli.T: threshold for posterior probability of toxicity for admissible set of toxicity;
## cutoff.eli.E: threshold for posterior probability of efficacy for admissible set of eficacy;
## ntrial: the number of simulated trial
## Uti: utility table for each possible outcome
## maxt: time assess windows for efficacy and toxicity 
## dist1:the underlying distribution of the time to toxicity/efficacy outcomes; dist1=1 is the uniform distribution, 
##        dist1=2 corresponds to the Weibull distribution.
## dist2: the underlying distribution of patient arrival time; dist3=1 is the uniform distribution, 
##        dist3=2 is the exponential distribution
## accrual: the accrual rate, i.e., the number of patients accrued in 1 unit of time
## alpha: a number from (0,1) that controls alpha*100% toxicity/efficacy events in (0, 1/2T)

#########################################################################################################



get.oc.UFO_B <- function(pT0.true,pT1.true,pE0.true,pE1.true,rho0,rho1,pI.true,
                        targetT=0.35,targetE=0.25,ncohort=20,cohortsize=3,
                        startdose=1,p.tox=0.358,cutoff.eli.T=0.95,
                        cutoff.eli.E=0.90,ntrial=5000,Uti,maxt=c(3,3),
                        dist1,dist2,accrual,alpha){
  library(Iso)
  # function to compute Pr(Tox=1,Eff=1)
  f1<-function(x,bn.m1,bn.m2,rho){
    ff<-dnorm(x,bn.m1,1)*(1-pnorm(0,bn.m2+rho*(x-bn.m1),sqrt(1-rho^2)))
    return(ff)
  }
  
  ndose=length(pI.true)
  
  joint.p=function(ttox,teff,c){
    ## ttox: marginal toxicity probability
    ## teff: marginal efficacy probability
    ## c: association parameter between tox and eff, c>0 indicates a positive correlation
    ndose=length(ttox) ## dose level
    out=matrix(rep(0,4*ndose),nrow=4) ## joint tox-eff probability matrix; row is dose level, column=(tox=0,eff=0; tox=0,eff=1; tox=1,eff=0; tox=1,eff=1)
    for (i in 1: ndose){
      out[4,i]=integrate(f1,bn.m1=qnorm(ttox[i]),bn.m2=qnorm(teff[i]),rho=c,lower=0,upper=Inf)$value
      out[2,i]=teff[i]-out[4,i]
      out[3,i]=ttox[i]-out[4,i]
      out[1,i]=1-out[4,i]-out[2,i]-out[3,i]
    }
    return(out)
  }
  
  pi0.j=joint.p(ttox=pT0.true, teff=pE0.true, c=rho0)
  pi1.j=joint.p(ttox=pT1.true, teff=pE1.true, c=rho1)
  
  ndose = length(pI.true)
  True_Uti = rep(0,ndose)
  for (tu in 1:ndose) {
    True_Uti[tu] = (1-pI.true[tu])*(pi0.j[1,tu]*Uti[1,1,1]+pi0.j[2,tu]*Uti[1,2,1]+
                                      pi0.j[3,tu]*Uti[2,1,1]+pi0.j[4,tu]*Uti[2,2,1])+
      pI.true[tu]*(pi1.j[1,tu]*Uti[1,1,2]+pi1.j[2,tu]*Uti[1,2,2]+pi1.j[3,tu]*Uti[2,1,2]+pi1.j[4,tu]*Uti[2,2,2])
    
  }
  #True_Uti
  
  
  ### data generation method: Gumbel copula
  Gumbel=function(ttox,teff,c){
    ## ttox: marginal toxicity probability
    ## teff: marginal efficacy probability
    ## c: association parameter between tox and eff, c>0 indicates a positive correlation
    ndose=length(ttox) ## dose level
    out=matrix(rep(0,4*ndose),nrow=4) ## joint tox-eff probability matrix; row is dose level, column=(tox=0,eff=0; tox=0,eff=1; tox=1,eff=0; tox=1,eff=1)
    for (j in 1: ndose){
      out[1,j]=(1-ttox[j])*(1-teff[j])+ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1 )/(exp(c)+1)
      out[2,j]=(1-ttox[j])*(teff[j])-ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1 )/(exp(c)+1)
      out[3,j]=(ttox[j])*(1-teff[j])-ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1 )/(exp(c)+1)
      out[4,j]=(ttox[j])*(teff[j])+ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1 )/(exp(c)+1)
    }
    return(out)
  }
  
  
  get.pts.outcome<-function(d,dist,cohortsize,pT0.true,pT1.true,pE0.true,pE1.true,timmu,pi0,pi1,maxt,alpha){
    
    d.pE0=pE0.true[d];d.pT0=pT0.true[d]
    d.pE1=pE1.true[d];d.pT1=pT1.true[d]
    
    pihalft_E0 = (1-alpha) * d.pE0
    pihalft_E1 = (1-alpha) * d.pE1
    pihalft_T0 = (1-alpha) * d.pT0
    pihalft_T1 = (1-alpha) * d.pT1
    
    library(truncdist)
    
    yI = rep(0,cohortsize)
    nT0 = rep(0,cohortsize)
    yT0 = rep(0,cohortsize)
    nT1 = rep(0,cohortsize)
    yT1 = rep(0,cohortsize)
    nE0 = rep(0,cohortsize)
    yE0 = rep(0,cohortsize)
    nE1 = rep(0,cohortsize)
    yE1 = rep(0,cohortsize)
    T0.time = rep(0,cohortsize)
    T1.time = rep(0,cohortsize)
    E0.time = rep(0,cohortsize)
    E1.time = rep(0,cohortsize)
    
    for (pts in 1:cohortsize) {
      yI[pts] = rbinom(1,1,timmu[d])
      if (yI[pts] == 0){
        res0=rmultinom(1,1,pi0[,d])
        yT0[pts] = res0[3]+res0[4]
        yE0[pts] = res0[2]+res0[4]
        nE0[pts] = 1 
        nT0[pts] = 1
        if(dist=="Weibull"){
          #Generate time to DLT
          T0.shape=log(log(1-d.pT0)/log(1-pihalft_T0))/log(2)
          T0.scale=maxt[1]/exp(log(-log(1-d.pT0))/T0.shape)
          #  print(paste("T.scale=",T.scale, "and T.shape=",T.shape))
          T0.time[pts]=ifelse(yT0[pts]==1,rtrunc(1, "weibull", a=0, b = maxt[1], shape=T0.shape,scale=T0.scale),maxt[1]+0.001);
          
          #Generate time to Efficacy
          E0.shape=log(log(1-d.pE0)/log(1-pihalft_E0))/log(2)
          E0.scale=maxt[2]/exp(log(-log(1-d.pE0))/E0.shape)
          #print(paste("E.scale=",E.scale, "and E.shape=",E.shape))
          E0.time[pts]=ifelse(yE0[pts]==1,rtrunc(1, "weibull", a=0, b = maxt[2], shape=E0.shape,scale=E0.scale),maxt[2]+0.001)
          
        }
        
        if(dist=="Uniform"){
          T0.time[pts]=ifelse(yT0[pts]==1, runif(1, 0, maxt[1]),maxt[1]+0.001)
          E0.time[pts]=ifelse(yE0[pts]==1, runif(1, 0, maxt[2]),maxt[2]+0.001)
        }
        
      }
      else{
        res1=rmultinom(1,1,pi1[,d])
        yT1[pts] = res1[3]+res1[4]
        yE1[pts] = res1[2]+res1[4]
        nE1[pts] = 1 
        nT1[pts] = 1
        
        if(dist=="Weibull"){
          #Generate time to DLT
          T1.shape=log(log(1-d.pT1)/log(1-pihalft_T1))/log(2)
          T1.scale=maxt[1]/exp(log(-log(1-d.pT1))/T1.shape)
          #  print(paste("T.scale=",T.scale, "and T.shape=",T.shape))
          T1.time[pts]=ifelse(yT1[pts]==1,rtrunc(1, "weibull", a=0, b = maxt[1], shape=T1.shape,scale=T1.scale),maxt[1]+0.001);
          
          #Generate time to Efficacy
          E1.shape=log(log(1-d.pE1)/log(1-pihalft_E1))/log(2)
          E1.scale=maxt[2]/exp(log(-log(1-d.pE1))/E1.shape)
          #print(paste("E.scale=",E.scale, "and E.shape=",E.shape))
          E1.time[pts]=ifelse(yE1[pts]==1,rtrunc(1, "weibull", a=0, b = maxt[2], shape=E1.shape,scale=E1.scale),maxt[2]+0.001)
          
        }
        
        if(dist=="Uniform"){
          T1.time[pts]=ifelse(yT1[pts]==1, runif(1, 0, maxt[1]),maxt[1]+0.001)
          E1.time[pts]=ifelse(yE1[pts]==1, runif(1, 0, maxt[2]),maxt[2]+0.001)
        }
      }
    }
    return(list(
      yI.i = yI,
      yI=sum(yI),
      nT=matrix(data = c(sum(nT0),sum(nT1)),nrow = 2),
      yT=matrix(data = c(sum(yT0),sum(yT1)),nrow = 2),
      nE=matrix(data = c(sum(nT0),sum(nE1)),nrow = 2),
      yE=matrix(data = c(sum(yE0),sum(yE1)),nrow = 2),
      T0.i = yT0[yI==0],
      T1.i = yT1[yI==1],
      E0.i = yE0[yI==0],
      E1.i = yE1[yI==1],
      T0.time = T0.time[T0.time!=0],
      E0.time = E0.time[E0.time!=0],
      T1.time = T1.time[T1.time!=0],
      E1.time = E1.time[E1.time!=0]
    ))
  }
  
  
  ####### Model for immune response
  model_I <- function(n, yi) {
    I_star = max(which(n!=0))
    I_star = ifelse(I_star==1,2,max(which(n!=0)))
    aic_I = rep(0, I_star)
    est_i=matrix(-0.5,nrow = I_star,ncol = I_star )
    for (l in 1:I_star) {
      if (l == 1) {
        sm = sum(yi[1:I_star])
        nm = sum(n[1:I_star])
        ql = sm / nm
        qfit = pava(y = ql, w = nm)
        ql_hat = qfit[1]
        est_i[l,1:I_star]=ql_hat
        likhood = (ql_hat ^ sm) * (1 - ql_hat) ^ (nm - sm)
        aic_I[l] = 2 * length(qfit) - 2 * log(likhood)
      } else{
        sk = (yi[1:I_star])[1:(l - 1)]
        nk = (n[1:I_star])[1:(l - 1)]
        qk = sk / nk
        qk[is.na(qk)] = 0
        sm = sum(yi[l:I_star])
        nm = sum(n[l:I_star])
        ql = sm / nm
        ql[is.na(ql)] = 0
        qfit = pava(y = c(qk, ql), w = c(nk, nm))
        qk_hat = qfit[1:(l - 1)]
        ql_hat = qfit[l]
        est_i[l,1:l]=c(qk_hat,ql_hat)
        est_i[l,l:I_star]=ql_hat
        likhood = prod((qk_hat ^ sk) * (1 - qk_hat) ^ (nk - sk)) * (ql_hat ^
                                                                      sm) * (1 - ql_hat) ^ (nm - sm)
        aic_I[l] = 2 * length(qfit) - 2 * log(likhood)
      }
    }
    list(AIC_I=aic_I,I_est=est_i)
  }
  
  
  ### Model for Toxicity
  model_T <- function(nti,yti,n) {
    T_star = max(which(n!=0))
    T_star = ifelse(T_star==1,2,max(which(n!=0)))
    nti = nti [,1:T_star]
    yti = yti [,1:T_star]
    pt = yti/nti
    pt[is.na(pt)] = 0
    T_est=biviso(y=pt,w=nti,fatal = F,warn = F)
    return(T_est)
  }
  ### Prepare data for Efficacy model
  modelE_dat <- function(nei,yei,n){
    E_star = max(which(n!=0))
    E_star = ifelse(E_star==1,2,max(which(n!=0)))
    nei = nei [,1:E_star]
    yei = yei [,1:E_star]
    pe1 = yei[1,]/nei[1,]
    nei2 = nei [2,]
    yei2 = yei [2,]
    pe2 = matrix(-10,nrow = E_star,ncol = E_star)
    wt_E = matrix(-10,nrow = E_star+1,ncol = E_star)
    wt_E [1,] = nei[1,]
    for (k in 1:E_star) {
      if (k == 1) {
        pe2 [k,] = c(rep(0,E_star-k),sum(yei2)/sum(nei2))
        wt_E [k+1,] = c(rep(0,E_star-k),sum(nei2))
      } else{
        ski = yei2[1:(k - 1)]
        nki = nei2[1:(k - 1)]
        qki = ski / nki
        smi = sum(yei2[k:E_star])
        nmi = sum(nei2[k:E_star])
        qli = smi / nmi
        pe2 [k,] = c(qki,rep(0,E_star-k),qli) 
        wt_E [k+1,] = c(nki,rep(0,E_star-k),nmi)
      }
    }
    pe1[is.na(pe1)] = 0
    pe2[is.na(pe2)] = 0
    return(list(est_E=rbind(pe1,pe2),wt_E=wt_E))
  }
  
  ### Estimate for Efficacy model
  modelE_est <- function(Est,Weight,n){
    E_star = max(which(n!=0))
    E_star = ifelse(E_star==1,2,max(which(n!=0)))
    Est_E = array(0,c(2,E_star,E_star))
    Weight_E = array(0,c(2,E_star,E_star))
    isoE = array(0,c(2,E_star,E_star))
    for (a in 1:E_star) {
      Est_E[1,,a] = Est[1,]
      Est_E[2,,a] = Est[a+1,]
      Weight_E[1,,a] = Weight[1,]
      Weight_E[2,,a] = Weight[a+1,]
    }
    for (b in 1:E_star) {
      if (b == E_star){
        isoE[,,b] = biviso(y=Est_E[,,b],w=Weight_E[,,b],fatal = F,warn = F)
      }
      else {
        isoE[,,b] = biviso(y=Est_E[,,b],w=Weight_E[,,b],fatal = F,warn = F)
        isoE[2,b:(E_star-1),b]= isoE[2,E_star,b]
      }
    }
    return(isoE)
  }
  
  
  ###AIC for Efficacy model
  modelE_AIC  <- function(nei,yei,n,pavaE) {
    E_star = max(which(n!=0))
    E_star = ifelse(E_star==1,2,max(which(n!=0)))
    aic_E = rep(0,E_star)
    for (c in 1:E_star) {
      v = as.vector((Est_E[,,c]^yei[,1:E_star])*(1-Est_E[,,c])^(nei[,1:E_star]-yei[,1:E_star]))
      aic_E[c] = 2*(E_star+c) -2*log(prod(v))
    }
    return(aic_E)
  }
  
  
  ###Final Estimate at each dose 
  Est_ETI = function(AICE,AICI,n,orderE,orderI,Est_T,Est_E,Est_I){
    jstar = max(which(n!=0))

    jstar = ifelse(jstar==1,2,max(which(n!=0)))
    
    if (jstar != 2) {
      weight=rep(0,6)
      Est_pai0 = matrix(0,nrow = 4,ncol = jstar)
      Est_pai1 = matrix(0,nrow = 4,ncol = jstar)
      
      AIC=c(AIC_I[orderI[1]]+AIC_E[orderE[1]],AIC_I[orderI[1]]+AIC_E[orderE[2]],
            AIC_I[orderI[2]]+AIC_E[orderE[1]],AIC_I[orderI[1]]+AIC_E[orderE[3]],
            AIC_I[orderI[2]]+AIC_E[orderE[2]],AIC_I[orderI[3]]+AIC_E[orderE[1]])
      weight=exp(-AIC/2)/sum(exp(-AIC/2))
      for (D in 1:jstar) {
        #pETI
        ##p000
        Est_pai0[1,D]=sum(c((1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[3]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[3],D])*(1-Est_E[1,D,orderE[1]])
        )*weight)*(1-Est_T[1,D])
        
        ##p001
        Est_pai1[1,D]=sum(c((Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[3]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[3],D])*(1-Est_E[2,D,orderE[1]])
        )*weight)*(1-Est_T[2,D])
        
        ##p010
        Est_pai0[2,D]=sum(c((1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[3]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[3],D])*(1-Est_E[1,D,orderE[1]])
        )*weight)*(Est_T[1,D])
        
        ##p011
        Est_pai1[2,D]=sum(c((Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[3]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[3],D])*(1-Est_E[2,D,orderE[1]])
        )*weight)*(Est_T[2,D])
        
        
        ##p100
        Est_pai0[3,D]=sum(c((1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[3]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[3],D])*(Est_E[1,D,orderE[1]])
        )*weight)*(1-Est_T[1,D])
        
        ##p101
        Est_pai1[3,D]=sum(c((Est_I[orderI[1],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[3]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[3],D])*(Est_E[2,D,orderE[1]])
        )*weight)*(1-Est_T[2,D])
        
        ##p110
        Est_pai0[4,D]=sum(c((1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[3]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[3],D])*(Est_E[1,D,orderE[1]])
        )*weight)*(Est_T[1,D])
        
        ##p111
        Est_pai1[4,D]=sum(c((Est_I[orderI[1],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[3]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[3],D])*(Est_E[2,D,orderE[1]])
        )*weight)*(Est_T[2,D])
        
      }
    } else{
      
      weight=rep(0,4)
      Est_pai0 = matrix(0,nrow = 4,ncol = jstar)
      Est_pai1 = matrix(0,nrow = 4,ncol = jstar)
      
      AIC=c(AIC_I[orderI[1]]+AIC_E[orderE[1]],AIC_I[orderI[1]]+AIC_E[orderE[2]],
            AIC_I[orderI[2]]+AIC_E[orderE[1]],AIC_I[orderI[2]]+AIC_E[orderE[2]])
      weight=exp(-AIC/2)/sum(exp(-AIC/2))
      for (D in 1:jstar) {
        #pETI
        ##p000
        Est_pai0[1,D]=sum(c((1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[2]])
        )*weight)*(1-Est_T[1,D])
        
        ##p001
        Est_pai1[1,D]=sum(c((Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[2]])
        )*weight)*(1-Est_T[2,D])
        
        ##p010
        Est_pai0[2,D]=sum(c((1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[2]])
        )*weight)*(Est_T[1,D])
        
        ##p011
        Est_pai1[2,D]=sum(c((Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[2]])
        )*weight)*(Est_T[2,D])
        
        
        ##p100
        Est_pai0[3,D]=sum(c((1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[2]])
        )*weight)*(1-Est_T[1,D])
        
        ##p101
        Est_pai1[3,D]=sum(c((Est_I[orderI[1],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[2]])
        )*weight)*(1-Est_T[2,D])
        
        ##p110
        Est_pai0[4,D]=sum(c((1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[2]])
        )*weight)*(Est_T[1,D])
        
        ##p111
        Est_pai1[4,D]=sum(c((Est_I[orderI[1],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[2]])
        )*weight)*(Est_T[2,D])
        
      }
      
    }
    return(list(weight=weight,Est_pai0=Est_pai0,Est_pai1=Est_pai1))
  }
  
  ###Est_Uti
  Est_Uti = function(Est0,Est1,Uti,n){
    jstar = max(which(n!=0))
    jstar = ifelse(jstar==1,2,max(which(n!=0)))
    Est_Uti=rep(0,jstar)
    for (et in 1:jstar) {
      Est_Uti[et] = Est0[1,et]*Uti[1,1,1]+Est0[2,et]*Uti[2,1,1]+
        Est0[3,et]*Uti[1,2,1]+Est0[4,et]*Uti[2,2,1]+
        Est1[1,et]*Uti[1,1,2]+Est1[2,et]*Uti[2,1,2]+
        Est1[3,et]*Uti[1,2,2]+Est1[4,et]*Uti[2,2,2]
    }
    return(Est_Uti=Est_Uti)
  }
  
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox,targetT,cutoff.eli.T) {
    at = 1 + ytox
    bt = 1 + n - ytox
    Tox_prob = 1 - pbeta(targetT, at, bt)
    if (length(Tox_prob[Tox_prob > cutoff.eli.T])>0){
      if (min(which(Tox_prob > cutoff.eli.T))==1){
        AT = 0
      } else{
        AT = c(1:(min(which(Tox_prob > cutoff.eli.T))-1))
      }} else {
        AT = c(1:length(n)) 
      }
    return(AT)
  }
  
  
  adm_eff <- function(n, yeff,targetE,cutoff.eli.E) {
    ae = 1 + yeff
    be = 1+ n - yeff
    Eff_prob = pbeta(targetE, ae, be)
    if (length(Eff_prob[Eff_prob > cutoff.eli.E])>0){
      AE = c((max(which(Eff_prob > cutoff.eli.E))+1):length(n))
    } else {
      AE = c(1:length(n)) 
    }
    return(AE)
  }
  
  
  
  dselect = rep(0, ntrial)
  ndose = length(pI.true)
  N = matrix(rep(100, ndose * ntrial), ncol = ndose)
  YTOX = matrix(rep(100, ndose * ntrial), ncol = ndose)
  YEFF = matrix(rep(100, ndose * ntrial), ncol = ndose)
  YI = matrix(rep(0, ndose * ntrial), ncol = ndose)
  pi0=Gumbel(ttox=pT0.true, teff=pE0.true, c=rho0)
  pi1=Gumbel(ttox=pT1.true, teff=pE1.true, c=rho1)
  
  durationV = rep(0, ntrial)
  
  for (trial in 1:ntrial) {
    # Immune response indicator for each subject
    yI.i = NULL
    #efficacy indicator for each subject
    yeff0.i = NULL  
    yeff1.i = NULL
    #toxicity indicator for each subject
    ytox0.i = NULL
    ytox1.i = NULL
    #dose for each subject
    dv0 = NULL
    dv1 = NULL
    dv = NULL
    
    # number of efficacy at each dose by interim time
    yI = rep(0,ndose)
    yeff0 = rep(0, ndose) 
    yeff1 = rep(0, ndose)
    n.eff0 = rep(0, ndose)    # Ess of efficacy at each dose
    n.eff1 = rep(0, ndose)
    yeff0.d = rep(0,ndose)
    yeff1.d = rep(0,ndose)
    n.eff0.d = rep(0,ndose)
    n.eff1.d = rep(0,ndose)
    
    
    ytox0 = rep(0, ndose) 
    ytox1 = rep(0, ndose)
    n.tox0 = rep(0, ndose)    # Ess of toxicity at each dose
    n.tox1 = rep(0, ndose)
    ytox0.d = rep(0,ndose)
    ytox1.d = rep(0,ndose)
    n.tox0.d = rep(0,ndose)
    n.tox1.d = rep(0,ndose)
    
    
    effobs = rep(0, ndose)
    toxobs = rep(0, ndose)
    
    ytox = rep(0, ndose)
    yeff = rep(0, ndose)
    n= rep(0, ndose)
    #n.d = rep(0,ndose)
    
    efft.enter = NULL
    efft0.event = NULL
    efft1.event = NULL
    efft.decision = 0
    
    toxt.enter = NULL
    toxt0.event = NULL
    toxt1.event = NULL
    toxt.decision = 0
    
    d = startdose ##starting dose level
    
    earlystop=0;              ## indicate whether the trial terminates early
    
    for (i in 1:ncohort) {
      if(d==0){break}
      # generate data for the new patient
      for (j in 1:cohortsize) {
        if (j == 1) {
          efft.enter = c(efft.enter, efft.decision)
          toxt.enter = efft.enter
        }
        else {
          if (dist2 == "Uniform") {
            efft.enter = c(efft.enter, efft.enter[length(efft.enter)] + runif(1, 0, 2 /
                                                                                accrual))
            toxt.enter = efft.enter
          }
          if (dist2 == "Exponential") {
            efft.enter = c(efft.enter, efft.enter[length(efft.enter)] + rexp(1, rate =
                                                                               accrual))
            toxt.enter = efft.enter
          }
        }
      }
      
      obscohort = get.pts.outcome(d=d,dist = dist1,cohortsize = cohortsize,pT0.true = pT0.true,
                                  pT1.true = pT1.true,pE0.true,pE1.true,timmu = pI.true,
                                  pi0 = pi0,pi1 = pi1,maxt=maxt,alpha = alpha)
      
      efft0.event = c(efft0.event, obscohort$E0.time)
      efft1.event = c(efft1.event, obscohort$E1.time)
      toxt0.event = c(toxt0.event, obscohort$T0.time)
      toxt1.event = c(toxt1.event, obscohort$T1.time)
      
      yI.i = c(yI.i, obscohort$yI.i)
      yeff0.i = c(yeff0.i, obscohort$E0.i)
      yeff1.i = c(yeff1.i, obscohort$E1.i)
      ytox0.i = c(ytox0.i, obscohort$T0.i)
      ytox1.i = c(ytox1.i, obscohort$T1.i)
      
      dv0 = c(dv0, rep(d, cohortsize-obscohort$yI))
      dv1 = c(dv1, rep(d, obscohort$yI))
      dv = c(dv, rep(d, cohortsize))
      
      
      efft.decision = efft.enter[length(efft.enter)]
      toxt.decision = toxt.enter[length(toxt.enter)]
      
      
      pending = 1
      
      while (pending == 1) {
        pending = 0
        if (i == ncohort) {
          efft.decision = efft.decision + maxt[2]
          toxt.decision = toxt.decision + maxt[1]
        }
        else {
          if (dist2 == "Uniform") {
            efft.decision = efft.decision + runif(1, 0, 2 / accrual)+0.001
            toxt.decision = efft.decision 
          }
          if (dist2 == "Exponential") {
            efft.decision = efft.decision + rexp(1, rate = accrual)+0.001
            toxt.decision = efft.decision
          }
        }
        
        efft0.enter = efft.enter[yI.i==0]
        efft1.enter = efft.enter[yI.i==1]
        
        toxt0.enter = toxt.enter[yI.i==0]
        toxt1.enter = toxt.enter[yI.i==1]
        
        # determine which observation are observed
        delta.eff0 = ((efft0.enter + efft0.event) <= efft.decision)
        delta.eff1 = ((efft1.enter + efft1.event) <= efft.decision)
        delta.tox0 = ((toxt0.enter + toxt0.event) <= toxt.decision)
        delta.tox1 = ((toxt1.enter + toxt1.event) <= toxt.decision)
        
        ## used for recording potential censoring time
        t.eff0 = pmin(efft0.event, efft.decision - efft0.enter, maxt[2])
        t.eff1 = pmin(efft1.event, efft.decision - efft1.enter, maxt[2])
        t.tox0 = pmin(toxt0.event, toxt.decision - toxt0.enter, maxt[1])
        t.tox1 = pmin(toxt1.event, toxt.decision - toxt1.enter, maxt[1])
        
        
        for(dd in 1:ndose){
          cset0 = (dv0 == dd)
          cset1 = (dv1 == dd)
          cset = (dv == dd)
          
          delta.eff0.curr = delta.eff0[cset0]
          delta.eff1.curr = delta.eff1[cset1]
          
          delta.tox0.curr = delta.tox0[cset0]
          delta.tox1.curr = delta.tox1[cset1]
          
          t.eff0.curr = t.eff0[cset0]
          t.eff1.curr = t.eff1[cset1]
          t.tox0.curr = t.tox0[cset0]
          t.tox1.curr = t.tox1[cset1]
          
          neff0.curr = sum((yeff0.i[cset0])[delta.eff0.curr == 1])
          neff1.curr = sum((yeff1.i[cset1])[delta.eff1.curr == 1])
          
          ntox0.curr = sum((ytox0.i[cset0])[delta.tox0.curr == 1])
          ntox1.curr = sum((ytox1.i[cset1])[delta.tox1.curr == 1])
          
          totalt.eff0 = t.eff0.curr[delta.eff0.curr == 0]
          totalt.eff1 = t.eff1.curr[delta.eff1.curr == 0]
          totalt.tox0 = t.tox0.curr[delta.tox0.curr == 0]
          totalt.tox1 = t.tox1.curr[delta.tox1.curr == 0]
          
          totalt.eff0 = sum(totalt.eff0) / maxt[2]
          totalt.eff1 = sum(totalt.eff1) / maxt[2]
          totalt.tox0 = sum(totalt.tox0) / maxt[1]
          totalt.tox1 = sum(totalt.tox1) / maxt[1]
          
          n.curr0 = sum(cset0)
          n.curr1 = sum(cset1)
          
          n.eff0.pend = sum(delta.eff0[cset0] == 0)
          n.eff1.pend = sum(delta.eff1[cset1] == 0)
          
          n.tox0.pend = sum(delta.tox0[cset0] == 0)
          n.tox1.pend = sum(delta.tox1[cset1] == 0)
          
          eff0obs = sum(delta.eff0[cset0])
          eff1obs = sum(delta.eff1[cset1])
          
          effobs[dd] = eff0obs + eff1obs
          
          tox0obs = sum(delta.tox0[cset0])
          tox1obs = sum(delta.tox1[cset1])
          
          toxobs[dd] = tox0obs + tox1obs
          
          n[dd] = sum(cset)
          
          n.eff0[dd] = n.curr0 - n.eff0.pend + totalt.eff0
          yeff0[dd] = neff0.curr
          n.eff1[dd] = n.curr1 - n.eff1.pend + totalt.eff1
          yeff1[dd] = neff1.curr
          
          n.tox0[dd] = n.curr0 - n.tox0.pend + totalt.tox0
          ytox0[dd] = ntox0.curr
          n.tox1[dd] = n.curr1 - n.tox1.pend + totalt.tox1
          ytox1[dd] = ntox1.curr
          
          yI[dd] = sum(yI.i[cset])
          
        }
        
        
        if ((effobs[d] > 0.5*n[d]) & (toxobs[d] >0.5* n[d])){
          pending=0
          ytox_I = matrix(rbind(ytox0,ytox1),nrow = 2)
          ytox = colSums(ytox_I)
          ntox_I =  matrix(rbind(n.tox0,n.tox1),nrow = 2)
          ntox = colSums(ntox_I)
          
          yeff_I = matrix(rbind(yeff0,yeff1),nrow = 2)
          yeff = colSums(yeff_I)
          neff_I =  matrix(rbind(n.eff0,n.eff1),nrow = 2)
          neff = colSums(neff_I)
          
          ### Get admissible set
          AT = adm_tox(n = ntox, ytox = ytox,targetT = targetT,cutoff.eli.T = cutoff.eli.T)
          AE = adm_eff(n = neff, yeff = yeff,targetE = targetE, cutoff.eli.E = cutoff.eli.E)
          
          A = intersect(AT,AE)
          
          if (length(A)==0){earlystop==1;d = 0;
          break}else{
            
            ## get observed toxcity rate
            pT_hat = ytox[d]/ntox[d]
            
            if (pT_hat>=p.tox & d==1){
              if (d %in% A){
                d_opt = d
              }else{
                earlystop = 1;d = 0;
                break}
            } else if (pT_hat>=p.tox & d!=1){
              if (length(A[A<=(d-1)])!=0){
                d_opt = max(A[A<=(d-1)])
              }else{
                if (d %in% A) {d_opt = d} else{earlystop = 1;d=0;break}
              }
            } else {
              if (d < ndose) {
                if (n[d+1]==0){
                  #d_opt = d+1
                  if ((d+1) %in% A){
                    d_opt = d+1} else {
                      d_opt = d
                    }
                } else{
                  admi_set=d;
                  if (d>1) {
                    if(length(A[A<=(d-1)])!=0){admi_set<-c(admi_set,max(A[A<=(d-1)]))} 
                  }
                  
                  if (length(A[A>=(d+1)])!=0){
                    admi_set<-c(admi_set,min(A[A>=(d+1)]))
                  }
                  admi_set = sort(admi_set)
                  
                  #######
                  model_I_result = model_I(n = n, yi = yI)
                  Est_T = model_T(nti = ntox_I, yti = ytox_I, n= n)
                  
                  modelE_data = modelE_dat (nei = neff_I,yei = yeff_I, n=n)
                  
                  Est_E = modelE_est(Est = modelE_data$est_E,
                                     Weight = modelE_data$wt_E,
                                     n=n)
                  
                  AIC_E=modelE_AIC(nei = neff_I,yei = yeff_I, n = n, pavaE = Est_E)
                  AIC_I=model_I_result$AIC_I
                  
                  orderE=order(AIC_E)
                  orderI=order(AIC_I)
                  
                  
                  Est_I=model_I_result$I_est
                  
                  Est_results = Est_ETI(AICE = AIC_E,AICI = AIC_I,n=n,orderE = orderE,
                                        orderI = orderI,Est_T = Est_T,Est_E = Est_E,
                                        Est_I = Est_I)
                  
                  Est_Utility = Est_Uti(Est0=Est_results$Est_pai0,Est1 = Est_results$Est_pai1,Uti = Uti,n=n)
                  #######
                  
                  d_opt = admi_set[which.max(Est_Utility[admi_set])]
                  
                }
              } else {
                admi_set=d;
                if (d>1) {
                  if(length(A[A<=(d-1)])!=0){admi_set<-c(admi_set,max(A[A<=(d-1)]))} 
                }
                admi_set = sort(admi_set)
                
                ######
                model_I_result = model_I(n = n, yi = yI)
                Est_T = model_T(nti = ntox_I, yti = ytox_I, n= n)
                
                modelE_data = modelE_dat (nei = neff_I,yei = yeff_I, n=n)
                
                Est_E = modelE_est(Est = modelE_data$est_E,
                                   Weight = modelE_data$wt_E,
                                   n=n)
                
                AIC_E=modelE_AIC(nei = neff_I,yei = yeff_I, n = n, pavaE = Est_E)
                AIC_I=model_I_result$AIC_I
                
                orderE=order(AIC_E)
                orderI=order(AIC_I)
                
                
                Est_I=model_I_result$I_est
                
                Est_results = Est_ETI(AICE = AIC_E,AICI = AIC_I,n=n,orderE = orderE,
                                      orderI = orderI,Est_T = Est_T,Est_E = Est_E,
                                      Est_I = Est_I)
                
                Est_Utility = Est_Uti(Est0=Est_results$Est_pai0,Est1 = Est_results$Est_pai1,Uti = Uti,n=n)
                ########
                
                d_opt = admi_set[which.max(Est_Utility[admi_set])]
              } 
            }
          }
          
          d<-d_opt
        }else{
          pending = 1;
        }
        
      }
    }
    
    YTOX[trial, ] = ytox
    YEFF[trial, ] = yeff
    YI[trial, ] = yI
    N[trial,]=n;
    durationV[trial] = max(toxt.decision,efft.decision)
    
    
    
    if (earlystop==0){ 
      pT = rep(targetT,length(n[n!=0]))
      pT = (ytox[n!=0])/(n[n!=0])
      pT = pava(y=pT,w=n[n!=0])
      d_mtd<-max(which(abs(pT-targetT)==min(abs(pT-targetT))))
      AT = adm_tox(n = n, ytox = ytox,targetT = targetT,cutoff.eli.T = cutoff.eli.T)
      AE = adm_eff(n = n, yeff = yeff,targetE = targetE, cutoff.eli.E = cutoff.eli.E)
      A = intersect(AT,AE)
      A = A[A<=d_mtd]
      
      if (length(A) == 0){
        dselect[trial] = 0
      } else {
        
        if (n[1]==60){
          dselect[trial] = 1
        } else{
          
          #######
          model_I_result = model_I(n = n, yi = yI)
          Est_T = model_T(nti = ntox_I, yti = ytox_I, n= n)
          
          modelE_data = modelE_dat (nei = neff_I,yei = yeff_I, n=n)
          
          Est_E = modelE_est(Est = modelE_data$est_E,
                             Weight = modelE_data$wt_E,
                             n=n)
          
          AIC_E=modelE_AIC(nei = neff_I,yei = yeff_I, n = n, pavaE = Est_E)
          AIC_I=model_I_result$AIC_I
          
          orderE=order(AIC_E)
          orderI=order(AIC_I)
          
          Est_I=model_I_result$I_est
          
          Est_results = Est_ETI(AICE = AIC_E,AICI = AIC_I,n=n,orderE = orderE,
                                orderI = orderI,Est_T = Est_T,Est_E = Est_E,
                                Est_I = Est_I)
          
          Est_Utility = Est_Uti(Est0=Est_results$Est_pai0,Est1 = Est_results$Est_pai1,Uti = Uti,n=n)
          #######
          dselect[trial] = A[which.max(Est_Utility[A])]
        }
      }
    } else{
      dselect[trial] = 0 
    }
  }
  selpercent = rep(0, ndose + 1)
  patpercent = matrix(rep(0, ntrial * ndose), ncol = ntrial, nrow = ndose)
  f <- function(x) {
    x[i] / sum(x)
  }
  ## Summarize results
  print("True Utility")
  cat(formatC(True_Uti, digits = 2, format = "f"),
      sep = " ", "\n")
  
  print("True Toxicity rate")
  cat(formatC(pT0.true*(1-pI.true)+pT1.true*(pI.true), digits = 2, format = "f"),
      sep = " ", "\n")
  
  print("True Efficacy rate")
  cat(formatC(pE0.true*(1-pI.true)+pE1.true*(pI.true), digits = 2, format = "f"),
      sep = " ", "\n")
  
  for (i in 0:ndose) {
    selpercent[(i + 1)] = sum(dselect == i) / ntrial * 100
  }
  print("selection probablity")
  cat(formatC(selpercent, digits = 1, format = "f"), sep = " ", "\n")
  for (i in 1:ndose) {
    patpercent[i, ] = apply(N, 1, f)
  }
  print("average percent of patients")
  cat(formatC(
    apply(patpercent, 1, mean) * 100,
    digits = 1,
    format = "f"
  ),
  sep = " ", "\n")
  print("average number of patients")
  cat(formatC(c(apply(N, 2, mean), sum(apply(
    N, 2, mean
  ))), digits = 1, format = "f"),
  sep = " ", "\n")
  
  print("average trial duration")
  cat(formatC(mean(durationV) , digits = 1, format = "f"),
      sep = " ", "\n")
  
}


#Utility
Uti <- array(0,c(2,2,2)) # order: tox, eff, immuno
Uti[,,1] <- matrix(c(0,0,80,35),nrow=2)
Uti[,,2] <- matrix(c(5,0,100,45),nrow=2)
Uti


##########################################################################################################
##########################################################################################################
## Function to get dose assignment for incoming cohort of the UFO design 
## To use this function, library "Iso" should be installed in R. 
##
##  Arguments:
## pts_data: patient-level data
## current_dose: current stayed dose
## total_dose: total number of under investigated doses
## targerT: highest acceptable toxicity rate
## targetE: lowest acceptable efficacy rate
## p.tox: dose deescalation boundary
## cutoff.eli.T: threshold for posterior probability of toxicity for admissible set of toxicity;
## cutoff.eli.E: threshold for posterior probability of efficacy for admissible set of eficacy;
## Uti: utility table for each possible outcome
## maxt: assessment window for efficacy/toxicity

get.df.UFO_B <- function(pts_data,current_dose, total_dose,
                       targetT=0.35,targetE=0.25,
                       p.tox=0.358,cutoff.eli.T=0.95,
                       cutoff.eli.E=0.90,Uti,maxt){
  library(Iso)
  ####### Model for immune response
  model_I <- function(n, yi) {
    I_star = max(which(n!=0))
    aic_I = rep(0, I_star)
    est_i=matrix(-0.5,nrow = I_star,ncol = I_star )
    for (l in 1:I_star) {
      if (l == 1) {
        sm = sum(yi[1:I_star])
        nm = sum(n[1:I_star])
        ql = sm / nm
        qfit = pava(y = ql, w = nm)
        ql_hat = qfit[1]
        est_i[l,1:I_star]=ql_hat
        likhood = (ql_hat ^ sm) * (1 - ql_hat) ^ (nm - sm)
        aic_I[l] = 2 * length(qfit) - 2 * log(likhood)
      } else{
        sk = (yi[1:I_star])[1:(l - 1)]
        nk = (n[1:I_star])[1:(l - 1)]
        qk = sk / nk
        qk[is.na(qk)] = 0
        sm = sum(yi[l:I_star])
        nm = sum(n[l:I_star])
        ql = sm / nm
        ql[is.na(ql)] = 0
        qfit = pava(y = c(qk, ql), w = c(nk, nm))
        qk_hat = qfit[1:(l - 1)]
        ql_hat = qfit[l]
        est_i[l,1:l]=c(qk_hat,ql_hat)
        est_i[l,l:I_star]=ql_hat
        likhood = prod((qk_hat ^ sk) * (1 - qk_hat) ^ (nk - sk)) * (ql_hat ^
                                                                      sm) * (1 - ql_hat) ^ (nm - sm)
        aic_I[l] = 2 * length(qfit) - 2 * log(likhood)
      }
    }
    list(AIC_I=aic_I,I_est=est_i)
  }
  
  
  ### Model for Toxicity
  model_T <- function(nti,yti,n) {
    T_star = max(which(n!=0))
    nti = nti [,1:T_star]
    yti = yti [,1:T_star]
    pt = yti/nti
    pt[is.na(pt)] = 0
    T_est=biviso(y=pt,w=nti,fatal = F,warn = F)
    return(T_est)
  }
  ### Prepare data for Efficacy model
  modelE_dat <- function(nei,yei,n){
    E_star = max(which(n!=0))
    nei = nei [,1:E_star]
    yei = yei [,1:E_star]
    pe1 = yei[1,]/nei[1,]
    nei2 = nei [2,]
    yei2 = yei [2,]
    pe2 = matrix(-10,nrow = E_star,ncol = E_star)
    wt_E = matrix(-10,nrow = E_star+1,ncol = E_star)
    wt_E [1,] = nei[1,]
    for (k in 1:E_star) {
      if (k == 1) {
        pe2 [k,] = c(rep(0,E_star-k),sum(yei2)/sum(nei2))
        wt_E [k+1,] = c(rep(0,E_star-k),sum(nei2))
      } else{
        ski = yei2[1:(k - 1)]
        nki = nei2[1:(k - 1)]
        qki = ski / nki
        smi = sum(yei2[k:E_star])
        nmi = sum(nei2[k:E_star])
        qli = smi / nmi
        pe2 [k,] = c(qki,rep(0,E_star-k),qli) 
        wt_E [k+1,] = c(nki,rep(0,E_star-k),nmi)
      }
    }
    pe1[is.na(pe1)] = 0
    pe2[is.na(pe2)] = 0
    return(list(est_E=rbind(pe1,pe2),wt_E=wt_E))
  }
  
  ### Estimate for Efficacy model
  modelE_est <- function(Est,Weight,n){
    E_star = max(which(n!=0))
    Est_E = array(0,c(2,E_star,E_star))
    Weight_E = array(0,c(2,E_star,E_star))
    isoE = array(0,c(2,E_star,E_star))
    for (a in 1:E_star) {
      Est_E[1,,a] = Est[1,]
      Est_E[2,,a] = Est[a+1,]
      Weight_E[1,,a] = Weight[1,]
      Weight_E[2,,a] = Weight[a+1,]
    }
    for (b in 1:E_star) {
      if (b == E_star){
        isoE[,,b] = biviso(y=Est_E[,,b],w=Weight_E[,,b],fatal = F,warn = F)
      }
      else {
        isoE[,,b] = biviso(y=Est_E[,,b],w=Weight_E[,,b],fatal = F,warn = F)
        isoE[2,b:(E_star-1),b]= isoE[2,E_star,b]
      }
    }
    return(isoE)
  }
  
  
  ###AIC for Efficacy model
  modelE_AIC  <- function(nei,yei,n,pavaE) {
    E_star = max(which(n!=0))
    aic_E = rep(0,E_star)
    for (c in 1:E_star) {
      v = as.vector((Est_E[,,c]^yei[,1:E_star])*(1-Est_E[,,c])^(nei[,1:E_star]-yei[,1:E_star]))
      aic_E[c] = 2*(E_star+c) -2*log(prod(v))
    }
    return(aic_E)
  }
  
  
  ###Final Estimate at each dose 
  Est_ETI = function(AICE,AICI,n,orderE,orderI,Est_T,Est_E,Est_I){
    jstar = max(which(n!=0))
    
    if (jstar != 2) {
      weight=rep(0,6)
      Est_pai0 = matrix(0,nrow = 4,ncol = jstar)
      Est_pai1 = matrix(0,nrow = 4,ncol = jstar)
      
      AIC=c(AIC_I[orderI[1]]+AIC_E[orderE[1]],AIC_I[orderI[1]]+AIC_E[orderE[2]],
            AIC_I[orderI[2]]+AIC_E[orderE[1]],AIC_I[orderI[1]]+AIC_E[orderE[3]],
            AIC_I[orderI[2]]+AIC_E[orderE[2]],AIC_I[orderI[3]]+AIC_E[orderE[1]])
      weight=exp(-AIC/2)/sum(exp(-AIC/2))
      for (D in 1:jstar) {
        #pETI
        ##p000
        Est_pai0[1,D]=sum(c((1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[3]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[3],D])*(1-Est_E[1,D,orderE[1]])
        )*weight)*(1-Est_T[1,D])
        
        ##p001
        Est_pai1[1,D]=sum(c((Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[3]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[3],D])*(1-Est_E[2,D,orderE[1]])
        )*weight)*(1-Est_T[2,D])
        
        ##p010
        Est_pai0[2,D]=sum(c((1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[3]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[3],D])*(1-Est_E[1,D,orderE[1]])
        )*weight)*(Est_T[1,D])
        
        ##p011
        Est_pai1[2,D]=sum(c((Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[3]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[3],D])*(1-Est_E[2,D,orderE[1]])
        )*weight)*(Est_T[2,D])
        
        
        ##p100
        Est_pai0[3,D]=sum(c((1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[3]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[3],D])*(Est_E[1,D,orderE[1]])
        )*weight)*(1-Est_T[1,D])
        
        ##p101
        Est_pai1[3,D]=sum(c((Est_I[orderI[1],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[3]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[3],D])*(Est_E[2,D,orderE[1]])
        )*weight)*(1-Est_T[2,D])
        
        ##p110
        Est_pai0[4,D]=sum(c((1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[3]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[3],D])*(Est_E[1,D,orderE[1]])
        )*weight)*(Est_T[1,D])
        
        ##p111
        Est_pai1[4,D]=sum(c((Est_I[orderI[1],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[3]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[3],D])*(Est_E[2,D,orderE[1]])
        )*weight)*(Est_T[2,D])
        
      }
    } else{
      
      weight=rep(0,4)
      Est_pai0 = matrix(0,nrow = 4,ncol = jstar)
      Est_pai1 = matrix(0,nrow = 4,ncol = jstar)
      
      AIC=c(AIC_I[orderI[1]]+AIC_E[orderE[1]],AIC_I[orderI[1]]+AIC_E[orderE[2]],
            AIC_I[orderI[2]]+AIC_E[orderE[1]],AIC_I[orderI[2]]+AIC_E[orderE[2]])
      weight=exp(-AIC/2)/sum(exp(-AIC/2))
      for (D in 1:jstar) {
        #pETI
        ##p000
        Est_pai0[1,D]=sum(c((1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[2]])
        )*weight)*(1-Est_T[1,D])
        
        ##p001
        Est_pai1[1,D]=sum(c((Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[2]])
        )*weight)*(1-Est_T[2,D])
        
        ##p010
        Est_pai0[2,D]=sum(c((1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(1-Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[2],D])*(1-Est_E[1,D,orderE[2]])
        )*weight)*(Est_T[1,D])
        
        ##p011
        Est_pai1[2,D]=sum(c((Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(1-Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[2],D])*(1-Est_E[2,D,orderE[2]])
        )*weight)*(Est_T[2,D])
        
        
        ##p100
        Est_pai0[3,D]=sum(c((1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[2]])
        )*weight)*(1-Est_T[1,D])
        
        ##p101
        Est_pai1[3,D]=sum(c((Est_I[orderI[1],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[2]])
        )*weight)*(1-Est_T[2,D])
        
        ##p110
        Est_pai0[4,D]=sum(c((1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[1],D])*(Est_E[1,D,orderE[2]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[1]]),
                            (1-Est_I[orderI[2],D])*(Est_E[1,D,orderE[2]])
        )*weight)*(Est_T[1,D])
        
        ##p111
        Est_pai1[4,D]=sum(c((Est_I[orderI[1],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[1],D])*(Est_E[2,D,orderE[2]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[1]]),
                            (Est_I[orderI[2],D])*(Est_E[2,D,orderE[2]])
        )*weight)*(Est_T[2,D])
        
      }
      
    }
    return(list(weight=weight,Est_pai0=Est_pai0,Est_pai1=Est_pai1))
  }
  
  ###Est_Uti
  Est_Uti = function(Est0,Est1,Uti,n){
    jstar = max(which(n!=0))
    Est_Uti=rep(0,jstar)
    for (et in 1:jstar) {
      Est_Uti[et] = Est0[1,et]*Uti[1,1,1]+Est0[2,et]*Uti[2,1,1]+
        Est0[3,et]*Uti[1,2,1]+Est0[4,et]*Uti[2,2,1]+
        Est1[1,et]*Uti[1,1,2]+Est1[2,et]*Uti[2,1,2]+
        Est1[3,et]*Uti[1,2,2]+Est1[4,et]*Uti[2,2,2]
    }
    return(Est_Uti=Est_Uti)
  }
  
  
  ###find admissible set of toxicity
  adm_tox <- function(n, ytox,targetT,cutoff.eli.T) {
    at = 1 + ytox
    bt = 1 + n - ytox
    Tox_prob = 1 - pbeta(targetT, at, bt)
    if (length(Tox_prob[Tox_prob > cutoff.eli.T])>0){
      if (min(which(Tox_prob > cutoff.eli.T))==1){
        AT = 0
      } else{
        AT = c(1:(min(which(Tox_prob > cutoff.eli.T))-1))
      }} else {
        AT = c(1:length(n)) 
      }
    return(AT)
  }
  
  
  adm_eff <- function(n, yeff,targetE,cutoff.eli.E) {
    ae = 1 + yeff
    be = 1+ n - yeff
    Eff_prob = pbeta(targetE, ae, be)
    if (length(Eff_prob[Eff_prob > cutoff.eli.E])>0){
      AE = c((max(which(Eff_prob > cutoff.eli.E))+1):length(n))
    } else {
      AE = c(1:length(n)) 
    }
    return(AE)
  }
  
  ####data####
  ndose = total_dose
  yI = rep(0,ndose)
  ytox = rep(0, ndose)
  ytox_I = matrix(0,ncol = ndose, nrow = 2)
  ntox_I = matrix(0,ncol = ndose, nrow = 2)
  yeff = rep(0, ndose)
  yeff_I = matrix(0,ncol = ndose, nrow = 2)
  neff_I = matrix(0,ncol = ndose, nrow = 2)
  ntox = rep(0, ndose)
  neff = rep(0, ndose)
  n = rep(0, ndose)
  earlystop=0; 
  
  
  output_function <- function(dat, maxt) {
    
    unique_doses <- unique(dat[, "dose"])
    n_doses <- length(unique_doses)
    
    # Matrix 1: Count of YT=1 for YI=0 and YI=1 at each dose
    matrix_YT <- matrix(0, 2, n_doses)
    
    for (i in 1:n_doses) {
      d <- unique_doses[i]
      matrix_YT[1, i] <- sum(dat[dat[, "YI"] == 0 & dat[, "dose"] == d, "YT"] == 1)
      matrix_YT[2, i] <- sum(dat[dat[, "YI"] == 1 & dat[, "dose"] == d, "YT"] == 1)
    }
    
    # Matrix 2: Count of YE=1 for YI=0 and YI=1 at each dose
    matrix_YE <- matrix(0, 2, n_doses)
    
    for (i in 1:n_doses) {
      d <- unique_doses[i]
      matrix_YE[1, i] <- sum(dat[dat[, "YI"] == 0 & dat[, "dose"] == d, "YE"] == 1)
      matrix_YE[2, i] <- sum(dat[dat[, "YI"] == 1 & dat[, "dose"] == d, "YE"] == 1)
    }
    
    # Vector: Count of YI=1 at each dose
    vector_YI <- integer(n_doses)
    
    for (i in 1:n_doses) {
      d <- unique_doses[i]
      vector_YI[i] <- sum(dat[dat[, "YI"] == 1 & dat[, "dose"] == d, "YI"])
    }
    
    
    # Vector1: Count of patients treated at each dose considering only YE
    vector_YE_counts <- numeric(n_doses)
    
    # Vector2: Count of patients treated at each dose considering only YT
    vector_YT_counts <- numeric(n_doses)
    
    # Matrices for counts of patients with YI=0 and YI=1 at each dose for YT and YE
    matrix_YI_YT <- matrix(0, 2, n_doses)
    matrix_YI_YE <- matrix(0, 2, n_doses)
    
    for (i in 1:n_doses) {
      d <- unique_doses[i]
      subrow=which(dat[,"dose"]==d)
      # For YE
      vector_YE_counts[i] <- sum(ifelse(dat[subrow, "YE"] == 1, 1, dat[subrow , "Follow-up time"] / maxt))
      
      
      # For YT
      vector_YT_counts[i] <- sum(ifelse(dat[subrow, "YT"] == 1, 1, dat[subrow , "Follow-up time"] / maxt))
      # For YI matrices
      for (j in 0:1) {  # loop through YI values
        #matrix_YI_YT[j+1, i] <- sum(ifelse(dat[subrow, "YI"] == j & dat[subrow, "YT"] == 1, 1,
        #                                   dat[intersect(subrow,subrow[which(dat[subrow, "YI"] ==j)]), "Follow-up time"] / maxt))
        
        matrix_YI_YT[j+1, i] <- sum(sum(dat[subrow, "YI"] == j & dat[subrow, "YT"]==1 ),
                                    ifelse(length(dat[subrow[dat[subrow, "YI"] == j & dat[subrow, "YT"] == 0],"Follow-up time"])==0,
                                           0, sum(dat[subrow[dat[subrow, "YI"] == j & dat[subrow, "YT"] == 0],"Follow-up time"]  / maxt))
                                           )
        
        
        matrix_YI_YE[j+1, i] <- sum(sum(dat[subrow, "YI"] == j & dat[subrow, "YE"] == 1),
                                     ifelse(length(dat[subrow[dat[subrow, "YI"] == j & dat[subrow, "YE"] == 0],"Follow-up time"])==0,
                                            0, sum(dat[subrow[dat[subrow, "YI"] == j & dat[subrow, "YE"] == 0],"Follow-up time"]  / maxt))
        )
        
      }
    }
    vector_dose_counts <- integer(n_doses)
    for (i in 1:n_doses) {
      d <- unique_doses[i]
      vector_dose_counts[i] <- sum(dat[, "dose"] == d)
    }
    return(list(matrix_YT = matrix_YT, matrix_YE = matrix_YE, vector_YI = vector_YI, 
                vector_YE_counts = vector_YE_counts, vector_YT_counts = vector_YT_counts, 
                matrix_YI_YT = matrix_YI_YT, matrix_YI_YE = matrix_YI_YE,
                vector_dose_counts = vector_dose_counts))
  }
  
  
  result <- output_function(dat = pts_data,maxt)
  result$matrix_YI_YT[is.na(result$matrix_YI_YT)]=0
  result$matrix_YI_YE[is.na(result$matrix_YI_YE)]=0
  
  high<-length(result$vector_dose_counts)
  
  yI[1:high] = result$vector_YI
  ytox[1:high] = colSums(result$matrix_YT)
  yeff[1:high] = colSums(result$matrix_YE)
  ntox[1:high] = result$vector_YT_counts
  neff[1:high] = result$vector_YE_counts
  ytox_I[,1:high] = result$matrix_YT
  yeff_I[,1:high] = result$matrix_YE
  ntox_I[,1:high] = result$matrix_YI_YT
  neff_I[,1:high] = result$matrix_YI_YE
  n[1:high] = result$vector_dose_counts
  
  
  ### Get admissible set
  AT = adm_tox(n = ntox, ytox = ytox,targetT = targetT,cutoff.eli.T = cutoff.eli.T)
  AE = adm_eff(n = neff, yeff = yeff,targetE = targetE, cutoff.eli.E = cutoff.eli.E)
  
  A = intersect(AT,AE)
  d =current_dose
  
  if (length(A)==0){
    d_opt = 0} else{
      ## get observed toxcity rate
      pT_hat = ytox[d]/n[d]
      if (pT_hat>=p.tox & d==1){
        if (d %in% A){
          d_opt = d
        }else{
          d_opt=0
        }
      } else if (pT_hat>=p.tox & d!=1){
        if (length(A[A<=(d-1)])!=0){
          d_opt = max(A[A<=(d-1)])
        }else{
          if (d %in% A) {d_opt = d} else{d_opt=0}
        }
      } else {
        if (d < ndose) {
          if (n[d+1]==0){
            #d_opt = d+1
            if ((d+1) %in% A){
              d_opt = d+1} else {
                d_opt = d
              }
          } else{
            admi_set=d;
            if (d>1) {
              if(length(A[A<=(d-1)])!=0){admi_set<-c(admi_set,max(A[A<=(d-1)]))} 
            }
            
            if (length(A[A>=(d+1)])!=0){
              admi_set<-c(admi_set,min(A[A>=(d+1)]))
            }
            admi_set = sort(admi_set)
            
            #######
            model_I_result = model_I(n = n, yi = yI)
            Est_T = model_T(nti = ntox_I, yti = ytox_I, n= n)
            
            modelE_data = modelE_dat (nei = neff_I,yei = yeff_I, n=n)
            
            Est_E = modelE_est(Est = modelE_data$est_E,
                               Weight = modelE_data$wt_E,
                               n=n)
            
            AIC_E=modelE_AIC(nei = neff_I,yei = yeff_I, n = n, pavaE = Est_E)
            AIC_I=model_I_result$AIC_I
            
            orderE=order(AIC_E)
            orderI=order(AIC_I)
            
            
            Est_I=model_I_result$I_est
            
            Est_results = Est_ETI(AICE = AIC_E,AICI = AIC_I,n=n,orderE = orderE,
                                  orderI = orderI,Est_T = Est_T,Est_E = Est_E,
                                  Est_I = Est_I)
            
            Est_Utility = Est_Uti(Est0=Est_results$Est_pai0,Est1 = Est_results$Est_pai1,Uti = Uti,n=n)
            #######
            
            d_opt = admi_set[which.max(Est_Utility[admi_set])]
          }
        } else {
          admi_set=d;
          if (d>1) {
            if(length(A[A<=(d-1)])!=0){admi_set<-c(admi_set,max(A[A<=(d-1)]))} 
          }
          admi_set = sort(admi_set)
          
          ######
          model_I_result = model_I(n = n, yi = yI)
          Est_T = model_T(nti = ntox_I, yti = ytox_I, n= n)
          
          modelE_data = modelE_dat (nei = neff_I,yei = yeff_I, n=n)
          
          Est_E = modelE_est(Est = modelE_data$est_E,
                             Weight = modelE_data$wt_E,
                             n=n)
          
          AIC_E=modelE_AIC(nei = neff_I,yei = yeff_I, n = n, pavaE = Est_E)
          AIC_I=model_I_result$AIC_I
          
          orderE=order(AIC_E)
          orderI=order(AIC_I)
          
          
          Est_I=model_I_result$I_est
          
          Est_results = Est_ETI(AICE = AIC_E,AICI = AIC_I,n=n,orderE = orderE,
                                orderI = orderI,Est_T = Est_T,Est_E = Est_E,
                                Est_I = Est_I)
          
          Est_Utility = Est_Uti(Est0=Est_results$Est_pai0,Est1 = Est_results$Est_pai1,Uti = Uti,n=n)
          ########
          
          d_opt = admi_set[which.max(Est_Utility[admi_set])]
        } 
      }
      
    }
  
  return(list("next dose"=d_opt))
  
}





########################################## Examples #######################################

####### generate operating characteristics of the UFO design ################################
# Assume that the target toxicity rate is 30%, the target efficacy rate is 25%, 
# and the sample size is 60 in cohorts of size 3.############################################
#### Scenario 1 ####
pI.true = c(0.7,0.7,0.7,0.7,0.7)
pT0.true = c(0.2,0.3,0.32,0.35,0.45)
pT1.true = c(0.22,0.4,0.55,0.6,0.6)
pE0.true = c(0.25,0.32,0.35,0.4,0.55)
pE1.true = c(0.6,0.6,0.6,0.6,0.6)
rho0=0
rho1=0

get.oc.UFO_B(pT0.true = pT0.true, pT1.true = pT1.true, pE0.true = pE0.true,
            pE1.true = pE1.true,pI.true = pI.true ,rho0 = rho0, rho1 = rho1, 
            Uti = Uti, ntrial = 10000,targetT = 0.3, targetE = 0.25,cutoff.eli.T = 0.9,
            cutoff.eli.E = 0.9, maxt=c(3,3),
            dist1="Weibull",dist2="Uniform",accrual=3,alpha = 0.5)
#[1] "True Utility"
#43.34 39.98 36.81 36.41 38.56 
#[1] "True Toxicity rate"
#0.21 0.37 0.48 0.52 0.56 
#[1] "True Efficacy rate"
#0.49 0.52 0.52 0.54 0.58 
#[1] "selection probablity"
#12.4 63.5 22.0 1.7 0.2 0.1 
#[1] "average percent of patients"
#57.6 27.9 9.7 3.8 1.0 
#[1] "average number of patients"
#31.0 14.6 4.6 1.7 0.6 52.4 
#[1] "average trial duration"
#30.5 

##### Scenario 4 #######
pI.true = c(0.35,0.6,0.6,0.6,0.6)
pT0.true = c(0.05,0.08,0.2,0.3,0.65)
pT1.true = c(0.06,0.14,0.34,0.6,0.7)
pE0.true = c(0.5,0.5,0.5,0.5,0.5)
pE1.true = c(0.6,0.6,0.6,0.6,0.6)
rho0=0
rho1=0

get.oc.UFO_B(pT0.true = pT0.true, pT1.true = pT1.true, pE0.true = pE0.true,
             pE1.true = pE1.true,pI.true = pI.true ,rho0 = rho0, rho1 = rho1, 
             Uti = Uti, ntrial = 10000,targetT = 0.3, targetE = 0.25,cutoff.eli.T = 0.9,
             cutoff.eli.E = 0.9, maxt=c(3,3),
             dist1="Weibull",dist2="Uniform",accrual=3,alpha = 0.5)
#[1] "True Utility"
#46.23 49.54 44.26 37.90 32.65 
#[1] "True Toxicity rate"
#0.05 0.12 0.28 0.48 0.68 
#[1] "True Efficacy rate"
#0.54 0.56 0.56 0.56 0.56 
#[1] "selection probablity"
#0.9 13.4 65.0 19.8 0.9 0.1 
#[1] "average percent of patients"
#16.5 49.4 24.5 7.5 2.1 
#[1] "average number of patients"
#9.4 29.6 14.7 4.5 1.3 59.5 
#[1] "average trial duration"
#37.7 


############  Hypothetical trial example from paper ########################################
############  Use get.df.UFO_B for dose assignment for incoming cohort #######################
############  Use get.dose_selection.UFO to select optimal dose by the end of trail ########
  
# Applying the function to the given data
### Cohort1 ###
dat <- matrix(c(0,0,0,0,0,1,0,1,0,90,90,80,1,1,1),3,5)
colnames(dat) <- c("YI","YT","YE","Follow-up time","dose")
get.df.UFO_B(pts_data=dat,current_dose=1, total_dose=5,
           targetT=0.3,targetE=0.25,
           p.tox=0.358,cutoff.eli.T=0.9,
           cutoff.eli.E=0.90,Uti=Uti,maxt = 90)
#$`next dose`
#[1] 2

### Cohort2 ###
dat <- matrix(c(0,0,0,1,0,1,0,0,1,0,0,0,0,1,0,0,1,0,
                90,90,90,90,90,80,1,1,1,2,2,2),6,5)
colnames(dat) <- c("YI","YT","YE","Follow-up time","dose")

get.df.UFO_B(pts_data=dat,current_dose=2, total_dose=5,
             targetT=0.3,targetE=0.25,
             p.tox=0.358,cutoff.eli.T=0.9,
             cutoff.eli.E=0.90,Uti=Uti,maxt = 90)
#$`next dose`
#[1] 3


### Cohort 4 ####
dat <- matrix(c(0,0,0,1,0,1,1,1,1,1,1,1,
                0,0,1,0,0,0,1,0,0,1,1,0,
                0,1,0,0,1,0,1,0,1,1,1,0,
                90,90,90,90,90,90,90,90,90,40,40,20,
                1,1,1,2,2,2,3,3,3,4,4,4),12,5)

colnames(dat) <- c("YI","YT","YE","Follow-up time","dose")

get.df.UFO_B(pts_data=dat,current_dose=4, total_dose=5,
             targetT=0.3,targetE=0.25,
             p.tox=0.358,cutoff.eli.T=0.9,
             cutoff.eli.E=0.90,Uti=Uti,maxt = 90)
#$`next dose`
#[1] 3

### Cohort 5 ####
dat <- matrix(c(0,0,0,1,0,1,1,1,1,1,1,1,1,0,1,
                0,0,1,0,0,0,1,0,0,1,1,1,0,1,0,
                0,1,0,0,1,0,1,0,1,1,1,1,0,1,0,
                90,90,90,90,90,90,90,90,90,90,90,90,
                90,40,70,
                1,1,1,2,2,2,3,3,3,4,4,4,3,3,3),15,5)

colnames(dat) <- c("YI","YT","YE","Follow-up time","dose")

get.df.UFO_B(pts_data=dat,current_dose=3, total_dose=5,
             targetT=0.3,targetE=0.25,
             p.tox=0.358,cutoff.eli.T=0.9,
             cutoff.eli.E=0.90,Uti=Uti,maxt = 90)
#$`next dose`
#[1] 2
