## get combined estimator by influence function
rm(list=ls())
setwd("/Users/chuanhong/Dropbox/00_Harvard/SurvCurvCalib/from_Chuan/code/to_share/")
library(locfit)
library(survival)
library(parallel)
source("Library.R")

rT.fun = function(nn){rexp(nn)}; rC.fun = function(nn){rexp(nn)}
get.semi<-function(Xi,Ci,Di,Zi,Ci.UL,Zi.UL){
  n = length(Ci); N = length(Ci.UL)
  ## ours
  Ni.t0  = sapply(1:n.t0,function(kk) I(Xi <= t0.all[kk])*Di)
  bet.t0 = sapply(1:n.t0,function(kk){
    ### Fit on labeled
    glm(Ni.t0[,kk]~Zi,family = binomial,weights = as.numeric(Ci>t0.all[kk]))$coef
  })
  
  WW.UL = cbind(1,Zi.UL)
  gg.UL  = sapply(1:n.t0,function(kk){
    tmp = expit(WW.UL%*%bet.t0[,kk])
    tmp
  })
  Semi = 1-sapply(1:n.t0,function(kk){
    sum(gg.UL[,kk]*{Ci.UL>t0.all[kk]})/sum(Ci.UL>t0.all[kk])
  })
  
  G_t  = sapply(t0.all,function(u) mean(Ci.UL>=u))
  tmp  = gg.UL*(1-gg.UL)
  tmp2 = tmp*outer(Ci.UL,t0.all,FUN=">")
  A    = sapply(1:n.t0,function(kk){
    tmp = t(WW.UL)%*%diag(tmp2[,kk])%*%WW.UL/N
    if(det(tmp)<1e-15) tmp = tmp+diag(rep(1e-15,2))
    solve(tmp)
  },simplify = "array")
  C    = t(tmp2)%*%WW.UL/N
  aa   = sapply(1:n.t0,function(kk) C[kk,]%*%A[,,kk]%*%t(cbind(1,Zi)*c({
    Ni.t0[,kk]-expit(cbind(1,Zi)%*%bet.t0[,kk])
  })*(Ci>t0.all[kk])))/VTM(G_t,n)
  
  return(list(Semi=Semi,Influence = aa))
}

n     = 200
N     = 2000
nboot = 500
lam   = 5
ncores  = 20
t0.all = seq(0.02,2,0.02); n.t0 = length(t0.all)
c.seq  = t0.all

set.seed(1234)
  Ti = rT.fun(n); ## failure time
  Ci = rC.fun(n); ## censoring time
  Zi = lam*log(Ti)+(1-lam)*log(Ci)+rnorm(n,sd=0.5);
  Xi = pmin(Ti,Ci); ## observed event time
  Di = I(Ti <= Ci) ## censoring indicator
  ## unlabeled
  Ti.UL = rT.fun(N); 
  Ci.UL = rC.fun(N);
  Zi.UL = lam*log(Ti.UL)+(1-lam)*log(Ci.UL)+rnorm(N,sd=0.5);
  Xi.UL = pmin(Ti.UL,Ci.UL); 
  Di.UL = I(Ti.UL <= Ci.UL)

  Semi = get.semi(Xi,Ci,Di,Zi,Ci.UL,Zi.UL)
 
