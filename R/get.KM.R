#' Semi-parametric estimator with influence function
#'
#' This function calculates the supervised KM estimator for survival data.
#'
#' @param Xi Vector of minimum of T and C for labeled data
#' @param Di Vector of event indicators for labeled data
#' @param t0.all Vector of time points at which the survival functions are estimated
#
#' @return A list containing the supervised KM estimator and the influence function
#' @export

get.KM<-function(Xi,Di, t0.all){
  Ni.t  = outer(Xi,t0.all,FUN = "<=")*Di
  wei   = Ni.t+outer(Xi,t0.all,FUN = ">")
  G_fit = survfit(Surv(Xi,1-Di)~1,type="kaplan-meier")
  G_hat = t(apply(outer(Xi,t0.all,FUN=pmin),1,
                  function(x) summary(G_fit,times=x,extend = TRUE)$surv))
  KM    = apply(wei/G_hat*outer(Xi,t0.all,FUN=">"),2,sum)/apply(wei/G_hat,2,sum)

  ## survival function of X
  ni  = sapply(Xi,function(u) sum(Xi>=u)) # survived at Xi
  SXX = ni/n # at each Xi (bar Y)
  ##
  ni  = 1/n/SXX^2*Di
  tmp = sapply(t0.all,function(u) outer(pmin(Xi,u),Xi,FUN=">=")%*%ni)
  Influence = Di/SXX*outer(Xi,t0.all,FUN="<=") - tmp
  Influence = Influence*VTM(KM,n)
  return(list(KM=KM,Influence=Influence))
}
