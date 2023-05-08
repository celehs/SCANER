
#' Non-parametric estimator
#'
#' This function calculates the non-parametric estimator for survival data
#'
#' @param Xi Vector of minimum of T and C for labeled data
#' @param Ci Vector of censoring times for labeled data
#' @param Di Vector of event indicators for labeled data
#' @param Zi Vector of covariates for labeled data
#' @param Ci.UL Vector of censoring times for unlabeled data
#' @param Zi.UL Vector of covariates for unlabeled data
#' @param h Bandwidth parameter for the kernel density estimation (default is set to 0.2)
#'
#' @return A list containing the non-parametric estimator
#' @export
get.NP<-function(Xi,Ci,Di,Zi,Ci.UL,Zi.UL,h=0.2){
  bet.t0   = sapply(1:n.t0,function(kk){
    ### Fit on labeled
    Ni.t0  = I(Xi <= t0.all[kk])*Di;
    glm(Ni.t0~Zi,family = binomial,weights = dnorm({log(Ci)-log(t0.all[kk])}/h)/h)$coef
  })

  tmp = cbind(1,Zi.UL)%*%bet.t0
  tmp = exp(tmp)/(1+exp(tmp))
  tmp2 = dnorm(outer(log(Ci.UL),log(t0.all),FUN="-")/h)/h
  NP = 1 - apply(tmp * tmp2,2,sum)/apply(tmp2,2,sum)

  return(NP)
}
