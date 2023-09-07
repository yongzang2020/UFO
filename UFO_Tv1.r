

## For targeted therapy


UFOT.OC=function(pitd,pie0t0d,pie0t1d,pie1t0d,pie1t1d,cohortsize,ncohort,phi.t,c.t,phi.e,c.e,omega,start,ntrial){
  ## pitd: vector of toxicity rate at different dose levels
  ## pie0t0d: vector for efficacy e=0 rate at different dose levels given t=0
  ## pie1t0d: vector for efficacy e=1 rate at different dose levels given t=0 
  ## pie0t1d: vector for efficacy e=0 rate at different dose levels given t=1
  ## pie1t1d: vector for efficacy e=1 rate at different dose levels given t=1
  ## omega: the utility value for (e=0,t=0;e=0,t=1;e=1,t=0;e=1,t=1;e=2,t=0;e=2,t=1)  
  ## efficacy is a categorical outcome: 2 for PR/CR; 1 for SD; 0 for DP.
  
  
  gen.P.T=function(pitd,pie0t0d,pie1t0d,pie0t1d,pie1t1d,omega){
    
    pie2t0d=1-pie0t0d-pie1t0d  
    pie2t1d=1-pie0t1d-pie1t1d
    
    
    ## marginal efficacy rates
    pie0d=pie0t0d*(1-pitd)+pie0t1d*pitd
    pie1d=pie1t0d*(1-pitd)+pie1t1d*pitd  
    pie2d=1-pie0d-pie1d 
    
    pie0t0=pie0t0d*(1-pitd) ## joint probability for e=0,t=0
    pie0t1=pie0t1d*pitd ## joint probability for e=0,t=1
    pie1t0=pie1t0d*(1-pitd) ## joint probability for e=1,t=0
    pie1t1=pie1t1d*pitd ## joint probability for e=1,t=1
    pie2t0=pie2t0d*(1-pitd) ## joint probability for e=2,t=0
    pie2t1=pie2t1d*pitd ## joint probability for e=2,t=1
    
    pi=rbind(pie0t0,pie0t1,pie1t0,pie1t1,pie2t0,pie2t1)
    ut=as.vector(t(as.matrix(omega))%*%pi) ## the utility at different dose levels
    return(list("pie0d"=pie0d, "pie1d"=pie1d, "pie2d"=pie2d, "pitd"=pitd, "pi"=pi, "ut"=ut    ))
    
    
  }
  
  
  ## estimate the joint probability of (e,t)
  est=function(n,omega){
   
    N=apply(n,2,sum) 
    ndose=dim(n)[2]
    pi=NULL
    for(i in 1:ndose){
      temp=n[,i]/N[i]
      pi=cbind(pi,temp)
    }
    
    
    ut=as.vector(t(as.matrix(omega))%*%pi) ## the utility at different dose levels
    return(list( "pi"=pi, "ut"=ut      ))
    
  }
  
  
  ## construct the admissible sets for toxicity and efficacy
  
  adm.tox=function(n,res,phi.t,c.t){
    p=NULL
    L=length(n)
    for (i in 1:L){
      p[i]=binom.test(res[i],n[i],phi.t,alternative="greater")[[3]]
    }
    if( any(p<c.t)==1  ){
      re=min( which(p<c.t)-1  )
    } else {re=L}
    return(re)
  }
  
  adm.eff=function(n,res,phi.e,c.e){
    p=NULL
    L=length(n)
    for (i in 1:L){
      p[i]=binom.test(res[i],n[i],phi.e, alternative="greater")[[3]]
    }
    re=which(p>=c.e)
    if( length(re)==0 ){
      return(0)
    } else {    return(re)   }
  }
  
  
  findobd=function(n,omega,admset){
    u=est(n,omega)$ut
    obd=which(u[admset]==max(u[admset]))[1]
    return(obd)
  }
  
  gen=gen.P.T(pitd,pie0t0d,pie1t0d,pie0t1d,pie1t1d,omega)
  t.pie0d=gen[[1]]
  t.pie1d=gen[[2]]
  t.pie2d=gen[[3]]
  t.pitd=gen[[4]]
  pi=gen[[5]]
  ut=gen[[6]]
  ndose=dim(pi)[2]
  
  N=matrix( rep(0,ndose*ntrial),ncol=ndose    )
  OBD=rep(0,ntrial)
  Toxrate=rep(0,ntrial)
  Eff0rate=rep(0,ntrial)
  Eff1rate=rep(0,ntrial)
  Eff2rate=rep(0,ntrial)
  
  for (i in 1: ntrial){
  
    y00=rep(0,ndose) ## e=0, t=0
    y01=rep(0,ndose) ## e=0, t=1
    y10=rep(0,ndose) ## e=1, t=0
    y11=rep(0,ndose) ## e=1, t=1
    y20=rep(0,ndose) ## e=2, t=0
    y21=rep(0,ndose) ## e=2, t=1
    n=rep(0,ndose)
    yt=rep(0,ndose)
    ye=rep(0,ndose)
    d=start
    estop=0
    for (j in 1: ncohort){
    
      n[d]=n[d]+cohortsize
      res=rmultinom(1,cohortsize,pi[,d])
      y00[d]=y00[d]+res[1]
      y01[d]=y01[d]+res[2]
      y10[d]=y10[d]+res[3]
      y11[d]=y11[d]+res[4]
      y20[d]=y20[d]+res[5]
      y21[d]=y21[d]+res[6]
      yt[d]=y01[d]+y11[d]+y21[d]
      ye[d]=y00[d]+y01[d]
      data=rbind(y00,y01,y10,y11,y20,y21)
      try=which(n>0) ## all the tired dose
      try.l=min(try) ## low bound for try
      try.u=max(try) ## up bound for try
      at=adm.tox(n[try],yt[try],phi.t,c.t)
      adm.high=at+try.l-1
      if (adm.high==0){
        estop=1
        break
      } else if (n[adm.high]==0){
        d=adm.high
      } else if ((adm.high==try.u)&(try.u<ndose)){
        d=adm.high+1
      } else {
        ae=adm.eff(n[try.l: adm.high], ye[try.l: adm.high], phi.e, c.e)
        if ((sum(ae)==0)&(try.l>1)){
          d=try.l-1
        } else if ((sum(ae)==0)&(try.l==1)){
          etsop=1
          break
        } else {
          temp=findobd(data[, try],omega,ae)
          d=ae[temp]+try.l-1
        }
      }
      
      
      
      
    }
    
    if (estop==1){
      dselect=0
    } else if (n[adm.high]==0){
      dselect=0
    } else if (sum(ae)==0){
      dselect=0
    } else{
      dselect=d
    }
    
    N[i,]=n
    OBD[i]=dselect
    Toxrate[i]=sum(yt)/sum(n)
    Eff0rate[i]=sum(ye)/sum(n)
    Eff1rate[i]=sum(y10+y11)/sum(n)
    Eff2rate[i]=sum(y20+y21)/sum(n)   
    
  }
  
  
  return(list("True Utility Value"=round(ut,digits=1),"True efficacy=0 rate"=t.pie0d,"True efficacy=1 rate"=t.pie1d, "True efficacy=2 rate"=t.pie2d, "True toxicity rate"=t.pitd,
              "OBD Selection Percentages"=round(table(OBD)*100/ntrial,digits=1) , "Average Number of Patients"=round(apply(N,2,mean),digits=1), "Average Total Number of Patients"=round(   sum( apply(N,2,mean)),digits=1),
              "Average Toxicity Rate"=round(mean(Toxrate)*100,digits=1), "Average disease progression Rate"=round(mean(Eff0rate)*100,digits=1), "Average stable disease Rate"=round(mean(Eff1rate)*100,digits=1), "Average PR/CR Rate"=round(mean(Eff2rate)*100,digits=1)  ))    
    
}  


pitd=c(0.15,0.15,0.15,0.4,0.5)
pie0t0d=c(0.5,0.4,0.3,0.2,0.1)
pie0t1d=c(0.4,0.3,0.2,0.1,0.05)
pie1t0d=c(0.3,0.2,0.1,0.05,0.05)
pie1t1d=c(0.2,0.1,0.05,0.05,0.05)
omega=c(10,0,50,30,100,80)
UFOT.OC(pitd,pie0t0d,pie0t1d,pie1t0d,pie1t1d,cohortsize=3,ncohort=20,phi.t=0.3,c.t=0.1,phi.e=0.5,c.e=0.1,omega,start=1,ntrial=1000)











