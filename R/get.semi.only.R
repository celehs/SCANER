## get our semiparametric estimator (without influence function variance)

## get our semiparametric estimator (without influence function variance)
get.semi.only<-function(Xi,Ci,Di,Zi,Ci.UL,Zi.UL, t0.all){
  n.t0=length(t0.all)
  ## ours
  Ni.t0  = sapply(1:n.t0,function(kk) I(Xi <= t0.all[kk])*Di)
  bet.t0 = sapply(1:n.t0,function(kk){
    ### Fit on labeled
    tryCatch(glm(Ni.t0[,kk]~Zi,family = binomial,weights = as.numeric(Ci>t0.all[kk]))$coef,error=function(e) NA)
  })
  bet.t0=do.call(cbind, bet.t0)

  WW.UL = cbind(1,Zi.UL)
  gg.UL  = sapply(1:n.t0,function(kk){
    if(any(is.na(bet.t0[,kk]))){tmp=rep(NA, dim(WW.UL)[1])}else{
      tmp = expit(WW.UL%*%bet.t0[,kk])}
    tmp
  })

  Semi = 1-sapply(1:n.t0,function(kk){
    sum(gg.UL[,kk]*{Ci.UL>t0.all[kk]})/sum(Ci.UL>t0.all[kk])
  })
  return(Semi)
}
