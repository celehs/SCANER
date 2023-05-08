#' Calculate AUC for the semi-parametric model
#'
#' This function calculates the area under the curve (AUC) for the
#' semi-parametric model using cross-validation.
#'
#' @param Xi A vector of event times for the labeled data.
#' @param Ci A vector of censoring times for the labeled data.
#' @param Di A vector of event indicators for the labeled data.
#' @param Zi A matrix of covariates for the labeled data.
#' @param Ci.UL A vector of censoring times for the unlabeled data.
#' @param Zi.UL A matrix of covariates for the unlabeled data.
#' @param t0.all A vector of time points at which to calculate the AUC.
#'
#' @return A vector containing the mean AUC for each time point in t0.all.
#' @export
get.semi.auc<-function(Xi,Ci,Di,Zi,Ci.UL,Zi.UL, t0.all){
  n.t0=length(t0.all)
  ## ours
  auc.t0=NULL
  for(iii in 1:10){
    id.v=sample(length(Xi), 0.5*length(Xi), replace=F)
    Xi.t=Xi[-id.v]
    Di.t=Di[-id.v]
    Ci.t=Ci[-id.v]
    Zi.t=Zi[-id.v,]
    Xi.v=Xi[id.v]
    Di.v=Di[id.v]
    Zi.v=Zi[id.v,]
    Ci.v=Ci[id.v]

    Ni.t0.t  = sapply(1:n.t0,function(kk) I(Xi.t <= t0.all[kk])*Di.t)
    Ni.t0.v  = sapply(1:n.t0,function(kk) I(Xi.v <= t0.all[kk])*Di.v)

    auc.t0 = rbind(auc.t0,sapply(1:n.t0,function(kk){
      ### Fit on labeled
      tryCatch({
        bb=glm(Ni.t0.t[,kk]~Zi.t,family = binomial,weights = as.numeric(Ci.t>t0.all[kk]))$coef
        ss=cbind(1,Zi.v)%*%bb
        ROC.Est.FUN(Ni.t0.v[,kk], ss, yy0=0.5, wgti=as.numeric(Ci.v>t0.all[kk]))[1]
      },error=function(e) NA)
    }))
  }
  auc.t=colMeans(auc.t0,na.rm=T)
  auc.t
}
