#' Compute Density Ratio estimator of survival function
#'
#' This function calculates the Density Ratio (DR) estimator of the survival function at the time points specified in `t0.all`.
#'
#' @param Zi Matrix of covariate values for the labeled data
#' @param Zi.UL Matrix of covariate values for the unlabeled data
#' @param Ti Vector of event times for the labeled data
#' @param Ci Vector of censoring times for the labeled data
#' @param Di Vector of the event indicators for the labeled data (1 for event, 0 for censored)
#' @param t0.all Vector of time points at which the survival function is to be estimated
#'
#' @return A list containing the Density Ratio estimator of the survival function (`DR`) at the specified time points in `t0.all`
#' @export

get.DR=function(Zi, Zi.UL, Ti, Ci, Di, t0.all){
  mygamma=Newton_DR(phi.t=Zi, phi.v=Zi.UL)
  DR=unlist(lapply(t0.all, function(tt) 1-sum(c(exp(t(mygamma)%*%t(Zi)))*(Ci>tt)*(Ti<tt))/sum(c(exp(t(mygamma)%*%t(Zi)))*(Ci>tt))))
  return(list(DR=DR))
}
