#' Generate synthetic data for survival analysis
#'
#' This function generates synthetic labeled and unlabeled data for survival analysis, simulating event times, censoring times, and covariates.
#'
#' @param n Integer specifying the number of labeled data points to generate
#' @param N Integer specifying the number of unlabeled data points to generate
#' @param lam Numeric value between 0 and 1 representing the proportion of event and censoring times used to generate the covariates
#'
#' @return A list containing the following generated data:
#' - `Ti`: Event times for the labeled data
#' - `Di`: Event indicators for the labeled data (1 for event, 0 for censored)
#' - `Ci`: Censoring times for the labeled data
#' - `Zi`: Covariate values for the labeled data
#' - `Xi`: Minimum of event and censoring times for the labeled data
#' - `Ti.UL`: Event times for the unlabeled data
#' - `Ci.UL`: Censoring times for the unlabeled data
#' - `Zi.UL`: Covariate values for the unlabeled data
#' - `Xi.UL`: Minimum of event and censoring times for the unlabeled data
#' @export

data_generation=function(n, N, lam){
Ti = rT.fun(n); Ci = rC.fun(n);
Zi = cbind(lam*log(Ti)+(1-lam)*log(Ci)+rnorm(n,sd=0.5),
           lam*log(Ti)+(1-lam)*log(Ci)+rnorm(n,sd=0.5))
Xi = pmin(Ti,Ci)
Di = I(Ti <= Ci)
## unlabeled
Ti.UL = rT.fun(N)
Ci.UL = rC.fun(N)
Zi.UL = cbind(lam*log(Ti.UL)+(1-lam)*log(Ci.UL)+rnorm(N,sd=0.5),
              lam*log(Ti.UL)+(1-lam)*log(Ci.UL)+rnorm(N,sd=0.5))
Xi.UL = pmin(Ti.UL,Ci.UL); Di.UL = I(Ti.UL <= Ci.UL)
return(list(Ti=Ti, Di=Di, Ci=Ci, Zi=Zi, Xi=Xi, Ti.UL=Ti.UL, Ci.UL=Ci.UL,Zi.UL=Zi.UL, Xi.UL=Xi.UL))
}
