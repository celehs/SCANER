#' Newton-DR Estimation for  Survival Data
#'
#' The `Newton_DR` function is an iterative algorithm to estimate the parameters in a
#' transformation model for survival data.
#'
#' @param phi.t A matrix with labeled data, where rows are observations and columns are covariates.
#' @param phi.v A matrix with unlabeled data, where rows are observations and columns are covariates.
#' @param weights An optional vector of weights for each observation. Default is 1 for equal weights.
#' @param max.iter Maximum number of iterations for the Newton-Raphson algorithm. Default is 100.
#' @param tol Tolerance for convergence of the Newton-Raphson algorithm. Default is 1e-5.
#' @param initial An optional vector of initial values for the algorithm. Default is a vector of zeros with length equal to the number of columns in phi.t.
#' @param lambda0 An optional initial value for the regularization parameter. If NULL, lambda0 is set to log(ncol(phi.t))/n.t^1.5. Default is NULL.
#'
#' @return A vector of parameter estimates for the transformation model.
#' @export

Newton_DR <- function(phi.t, phi.v, weights = 1,
                      max.iter = 100, tol = 1e-5, initial = rep(0, ncol(phi.t)),
                      lambda0 = NULL){

  n.t <- nrow(phi.t);
  n.v <- nrow(phi.v);
  if(is.null(lambda0)){lambda0 = log(ncol(phi.t))/n.t^1.5};

  error <- Inf
  iter <- 0
  gamma <- initial
  phi.v.mean <- rowMeans(t(phi.v))

  while(iter < max.iter & error > tol){
    gamma_old <- gamma
    z <- as.vector(phi.t %*% gamma)
    w <- exp(z) * weights
    phiT.phi <- crossprod(phi.t, w * phi.t) / n.t
    E.phi <- t(phi.t) %*% w / n.t

    gamma <- gamma + solve(phiT.phi + diag(rep(lambda0, ncol(phiT.phi)))) %*% (phi.v.mean - E.phi)
    error <- sqrt(mean((gamma - gamma_old)^2))
    #print(error)
  }
  return(gamma)

}
