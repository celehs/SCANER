#' Calculate Deep Learning-based estimator
#'
#' This function computes the deep learning-based estimator for survival probabilities
#' at each specified time point using a deep survival model.
#'
#' @param Zi A matrix of covariates for the labeled data
#' @param Zi.UL A matrix of covariates for the unlabeled data
#' @param Xi A numeric vector representing event/censoring times for the labeled data
#' @param Ci A numeric vector representing censoring times for the labeled data
#' @param Di A numeric vector representing event indicators for the labeled data
#' @param t0.all A numeric vector representing the time points at which to estimate survival probabilities
#' @return A data frame containing the time points and the corresponding deep learning-based estimates of survival probabilities
#' @export

get.DL=function(Zi, Zi.UL, Xi, Ci, Di, t0.all){
  mydat=data.frame(time=Xi, status=Di, Zi)
  mydat2=data.frame(time=Xi, status=Di, Zi)

  fit=deepsurv(data = mydat, frac = 0.3, activation = "relu",
               num_nodes = c(4L, 8L, 4L, 2L), dropout = 0.1, early_stopping = TRUE, epochs = 100L,
               batch_size = 32L)
  colnames(Zi.UL)=colnames(Zi)
  p.pat=predict(fit, data.frame(Zi.UL))
  t0.all.new=colnames(p.pat)
  n.t0.new=length(t0.all.new)
  pp = sapply(1:n.t0.new,function(kk){
    sum(p.pat[,kk]*{Ci.UL>t0.all.new[kk]})/sum(Ci.UL>t0.all.new[kk])
  })
  DL=colMeans(p.pat)
  DL=data.frame(t0=names(DL), DL)
  row.names(DL)=NULL
  return(DL)
}
