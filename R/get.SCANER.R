#' Calculate SCANER estimator
#'
#' This function computes the combined SCANER estimator for survival probabilities
#' at each specified time point using both the semi-supervised and supervised estimators.
#'
#' @param n.t0 Integer representing the number of time points to evaluate
#' @param n Integer representing the sample size
#' @param Semi A list containing the semi-supervised estimates and their influence values
#' @param KM A list containing the Kaplan-Meier estimates and their influence values
#' @return A numeric vector of SCANER estimates for survival probabilities at each time point
#' @export

get.SCANER=function(n.t0, n, Semi, KM){
Comb.Cov = array(0,dim=c(2,2,n.t0))
Comb.Cov[1,1,] = apply(Semi$Influence,2,var)/n
Comb.Cov[2,2,] = apply(KM$Influence,2,var)/n
Comb.Cov[1,2,] = Comb.Cov[2,1,] = sapply(1:n.t0,function(kk) cov(Semi$Influence[,kk],KM$Influence[,kk]))/n

## inverse
Comb.Cov = sapply(1:n.t0,function(kk){
  if(any(is.na(Comb.Cov[,,kk]))) return(matrix(NA,ncol=2,nrow=2))
  if(det(Comb.Cov[,,kk])<1e-15) Comb.Cov[,,kk]=Comb.Cov[,,kk]+diag(rep(1e-15,2))
  return(solve(Comb.Cov[,,kk]))
},simplify = "array")

## combined estimator
Comb = sapply(1:n.t0,function(kk) t(c(1,1))%*%Comb.Cov[,,kk]%*%c(Semi$Semi[kk],KM$KM[kk])/sum(Comb.Cov[,,kk]))
Comb
}
