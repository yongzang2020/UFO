## For Immunotherapy 






UFOI.OC=function(piid,piei,pit0d,pit1d,cohortsize,ncohort,phi.t,c.t,phi.e,c.e,omega,start,ntrial){
  ## piid: vector of immune response rate at different dose levels
  ## piei: matrix for efficacy rate given immune rate: rows for e=0,1,2; columns for i=0,1.
  ## pit0d: vector for toxicity rate at different dose levels given i=0
  ## pit1d: vector for toxicity rate at different dose levels given i=1 
  ## efficacy is a categorical outcome: 2 for PR/CR; 1 for SD; 0 for DP.
  
  
  ## determine the joint toxicity, efficacy and immune response rates 
  gen.P.I=function(piid,piei,pit0d,pit1d){
    
    pie0d=piei[1,1]*(1-piid)+piei[1,2]*piid ## probability of e=0 at different dose levels
    pie1d=piei[2,1]*(1-piid)+piei[2,2]*piid ## probability of e=1 at different dose levels
    pie2d=piei[3,1]*(1-piid)+piei[3,2]*piid ## probability of e=2 at different dose levels  
    pitd=(1-piid)*pit0d+piid*pit1d## toxicity rates at different dose levels 
    
    
    pii0e0t0=(1-piid)*piei[1,1]*(1-pit0d)
    pii0e0t1=(1-piid)*piei[1,1]*pit0d
    pii0e1t0=(1-piid)*piei[2,1]*(1-pit0d)
    pii0e1t1=(1-piid)*piei[2,1]*pit0d
    pii0e2t0=(1-piid)*piei[3,1]*(1-pit0d)
    pii0e2t1=(1-piid)*piei[3,1]*pit0d
    
    pii1e0t0=(piid)*piei[1,2]*(1-pit1d)
    pii1e0t1=(piid)*piei[1,2]*pit1d
    pii1e1t0=(piid)*piei[2,2]*(1-pit1d)
    pii1e1t1=(piid)*piei[2,2]*pit1d
    pii1e2t0=(piid)*piei[3,2]*(1-pit1d)
    pii1e2t1=(piid)*piei[3,2]*pit1d
    
    pi=rbind(pii0e0t0,pii0e0t1,pii0e1t0,pii0e1t1,pii0e2t0,pii0e2t1, pii1e0t0,pii1e0t1,pii1e1t0,pii1e1t1,pii1e2t0,pii1e2t1  ) ## joint probability of i, e, t at different dose levels
    return(list("piid"=piid, "pie0d"=pie0d, "pie1d"=pie1d, "pie2d"=pie2d, "pitd"=pitd, "pi"=pi    ))
    
  }
  
  
  
  
  ## estimate the joint probability of (i,e,t)
  est=function(n,omega){
    piid=apply(n[7:12,],2,sum)/apply(n,2,sum) ## estimator for probability of i=1
    pit0d=NULL
    pit1d=NULL
    for (i in 1: length(piid)){
      if (  ( piid[i]==0   )|(piid[i]==1)    ){
        pit0d[i]=pit1d[i]=( n[2,i]+n[4,i]+n[6,i]+n[8,i]+n[10,i]+n[12,i]    )/( sum(n[,i])     )
      }else{
        pit0d[i]=( n[2,i]+n[4,i]+n[6,i]    )/( sum(n[1:6,i])     )
        pit1d[i]=( n[8,i]+n[10,i]+n[12,i]    )/( sum(n[7:12,i])     ) 
      }
    }
    n.comb=apply(n,1,sum) ## combine dose levels
    pii=sum(n.comb[7:12])/sum(n.comb)
    pie0i0=NULL
    pie0i1=NULL
    pie1i0=NULL
    pie1i1=NULL
    pie2i0=NULL
    pie2i1=NULL
    if( ( pii==0   )|(  pii==1  )      ){
      pie0i0=pie0i1=( n.comb[1]+n.comb[2]+n.comb[7]+n.comb[8])/sum(n.comb)
      pie1i0=pie1i1=( n.comb[3]+n.comb[4]+n.comb[9]+n.comb[10])/sum(n.comb)
      pie2i0=pie2i1=( n.comb[5]+n.comb[6]+n.comb[11]+n.comb[12])/sum(n.comb)
    }else{
      pie0i0=( n.comb[1]+n.comb[2])/sum(n.comb[1:6])
      pie0i1=( n.comb[7]+n.comb[8])/sum(n.comb[7:12])
      pie1i0=( n.comb[3]+n.comb[4])/sum(n.comb[1:6])
      pie1i1=( n.comb[9]+n.comb[10])/sum(n.comb[7:12])
      pie2i0=( n.comb[5]+n.comb[6])/sum(n.comb[1:6])
      pie2i1=( n.comb[11]+n.comb[12])/sum(n.comb[7:12])    
    }
    
    
    pii0e0t0=(1-piid)*pie0i0*(1-pit0d)
    pii0e0t1=(1-piid)*pie0i0*pit0d
    pii0e1t0=(1-piid)*pie1i0*(1-pit0d)
    pii0e1t1=(1-piid)*pie1i0*pit0d
    pii0e2t0=(1-piid)*pie2i0*(1-pit0d)
    pii0e2t1=(1-piid)*pie2i0*pit0d
    
    pii1e0t0=(piid)*pie0i1*(1-pit1d)
    pii1e0t1=(piid)*pie0i1*pit1d
    pii1e1t0=(piid)*pie1i1*(1-pit1d)
    pii1e1t1=(piid)*pie1i1*pit1d
    pii1e2t0=(piid)*pie2i1*(1-pit1d)
    pii1e2t1=(piid)*pie2i1*pit1d
    
    pi=rbind(pii0e0t0,pii0e0t1,pii0e1t0,pii0e1t1,pii0e2t0,pii0e2t1, pii1e0t0,pii1e0t1,pii1e1t0,pii1e1t1,pii1e2t0,pii1e2t1  ) 
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
  
  
  gen=gen.P.I(piid,piei,pit0d,pit1d)
  t.piid=gen[[1]]
  t.pie0d=gen[[2]]
  t.pie1d=gen[[3]]
  t.pie2d=gen[[4]]
  t.pitd=gen[[5]]
  pi=gen[[6]]
  
  ndose=dim(pi)[2]
  
  
  tu=NULL
  for(i in 1: ndose){
    tu[i]=sum(pi[,i]*omega)
  }
  
  N=matrix( rep(0,ndose*ntrial),ncol=ndose    )
  OBD=rep(0,ntrial)
  Toxrate=rep(0,ntrial)
  Eff0rate=rep(0,ntrial)
  Eff1rate=rep(0,ntrial)
  Eff2rate=rep(0,ntrial)
  Imurate=rep(0,ntrial)
  
  for (i in 1: ntrial){
    
    y000=rep(0,ndose) ## i=0, e=0, t=0
    y001=rep(0,ndose) ## i=0, e=0, t=1
    y010=rep(0,ndose) ## i=0, e=1, t=0
    y011=rep(0,ndose) ## i=0, e=1, t=1
    y020=rep(0,ndose) ## i=0, e=2, t=0
    y021=rep(0,ndose) ## i=0, e=2, t=1
    y100=rep(0,ndose) ## i=1, e=0, t=0
    y101=rep(0,ndose) ## i=1, e=0, t=1
    y110=rep(0,ndose) ## i=1, e=1, t=0
    y111=rep(0,ndose) ## i=1, e=1, t=1
    y120=rep(0,ndose) ## i=1, e=2, t=0
    y121=rep(0,ndose) ## i=1, e=2, t=1
    n=rep(0,ndose)
    yt=rep(0,ndose)
    ye=rep(0,ndose)
    d=start
    estop=0
    for (j in 1: ncohort){
      n[d]=n[d]+cohortsize
      res=rmultinom(1,cohortsize,pi[,d])
      y000[d]=y000[d]+res[1]
      y001[d]=y001[d]+res[2]
      y010[d]=y010[d]+res[3]
      y011[d]=y011[d]+res[4]
      y020[d]=y020[d]+res[5]
      y021[d]=y021[d]+res[6]
      y100[d]=y100[d]+res[7]
      y101[d]=y101[d]+res[8]
      y110[d]=y110[d]+res[9]
      y111[d]=y111[d]+res[10]
      y120[d]=y120[d]+res[11]
      y121[d]=y121[d]+res[12]
      yt[d]=y001[d]+y011[d]+y021[d]+y101[d]+y111[d]+y121[d]
      ye[d]=y000[d]+y001[d]+y100[d]+y101[d]
      data=rbind(y000,y001,y010,y011,y020,y021,y100,y101,y110,y111,y120,y121)
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
    Eff1rate[i]=sum(y010+y011+y110+y111)/sum(n)
    Eff2rate[i]=sum(y020+y021+y120+y121)/sum(n)
    Imurate[i]=sum(y100+y101+y110+y111+y120+y121)/sum(n)
  }
  
  return(list("True Utility Value"=round(tu,digits=1),"True immuno response rate"=t.piid,"True efficacy=0 rate"=t.pie0d,"True efficacy=1 rate"=t.pie1d, "True efficacy=2 rate"=t.pie2d, "True toxicity rate"=t.pitd,
              "OBD Selection Percentages"=round(table(OBD)*100/ntrial,digits=1) , "Average Number of Patients"=round(apply(N,2,mean),digits=1), "Average Total Number of Patients"=round(   sum( apply(N,2,mean)),digits=1),
              "Average Toxicity Rate"=round(mean(Toxrate)*100,digits=1), "Average disease progression Rate"=round(mean(Eff0rate)*100,digits=1), "Average stable disease Rate"=round(mean(Eff1rate)*100,digits=1), "Average PR/CR Rate"=round(mean(Eff2rate)*100,digits=1), "Average Immuno Rate"=round(mean(Imurate)*100,digits=1)  ))    


}


# piid=c(0.2,0.3,0.4,0.5,0.6)
# piei=matrix( c(0.8,0.1,0.1,0.9,0.05,0.05), nrow=3   )
# pit0d=seq(0.05,0.25,0.05)
# pit1d=seq(0.1,0.3,0.05)
# omega=c(0, 0, 50,10,80,35,5,0,70,20,100,45)
# UFOI.OC(piid,piei,pit0d,pit1d,cohortsize=3,ncohort=20,phi.t=0.3,c.t=0.1,phi.e=0.5,c.e=0.1,omega,start=1,ntrial=1000)


# piid=c(0.1,0.2,0.3,0.4,0.8)
# piei=matrix( c(0.2,0.5,0.3,0.1,0.45,0.45), nrow=3   )
# pit0d=c(0.5,0.6,0.7,0.8,0.9)
# pit1d=c(0.5,0.6,0.7,0.8,0.9)
# omega=c(0, 0, 50,10,80,35,5,0,70,20,100,45)
# UFOI.OC(piid,piei,pit0d,pit1d,cohortsize=3,ncohort=20,phi.t=0.3,c.t=0.1,phi.e=0.5,c.e=0.1,omega,start=1,ntrial=1000)
# # 
# 
# piid=c(0.2,0.3,0.4,0.5,0.6)
# piei=matrix( c(0.1,0.5,0.4,0.2,0.4,0.4), nrow=3   )
# pit0d=seq(0.25,0.05,-0.05)
# pit1d=seq(0.3,0.1,-0.05)
# omega=c(0, 0, 50,10,80,35,5,0,70,20,100,45)
# UFOI.OC(piid,piei,pit0d,pit1d,cohortsize=3,ncohort=20,phi.t=0.3,c.t=0.1,phi.e=0.5,c.e=0.1,omega,start=1,ntrial=1000)
# 


