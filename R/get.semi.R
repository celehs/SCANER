
#' Semi-supervised estimator with influence function
#'
#' This function calculates the semi-supervised estimator for survival data.
#'
#' @param Xi Vector of minimum of T and C for labeled data
#' @param Ci Vector of censoring times for labeled data
#' @param Di Vector of event indicators for labeled data
#' @param Zi Vector of covariates for labeled data
#' @param Ci.UL Vector of censoring times for unlabeled data
#' @param Zi.UL Vector of covariates for unlabeled data
#'
#' @return A list containing the semi-parametric estimator and the influence function
#' @export

get.semi<-function(Xi,Ci,Di,Zi,Ci.UL,Zi.UL, t0.all){
  n = length(Ci); N = length(Ci.UL)
  ## ours
  Ni.t0  = sapply(1:n.t0,function(kk) I(Xi <= t0.all[kk])*Di)
  bet.t0 = sapply(1:n.t0,function(kk){
    ### Fit on labeled
    tryCatch(glm(Ni.t0[,kk]~Zi,family = binomial,weights = as.numeric(Ci>t0.all[kk]))$coef, error=function(e) NA)
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
    if(det(tmp)<1e-15) tmp = tmp+diag(rep(1e-15,dim(tmp)[2]))
    solve(tmp)
  },simplify = "array")
  C    = t(tmp2)%*%WW.UL/N
  aa   = sapply(1:n.t0,function(kk) C[kk,]%*%A[,,kk]%*%t(cbind(1,Zi)*c({
    Ni.t0[,kk]-expit(cbind(1,Zi)%*%bet.t0[,kk])
  })*(Ci>t0.all[kk])))/VTM(G_t,n)

  return(list(Semi=Semi,Influence = aa))
}

