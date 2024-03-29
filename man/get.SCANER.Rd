% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get.SCANER.R
\name{get.SCANER}
\alias{get.SCANER}
\title{Calculate SCANER estimator}
\usage{
get.SCANER(n.t0, n, Semi, KM)
}
\arguments{
\item{n.t0}{Integer representing the number of time points to evaluate}

\item{n}{Integer representing the sample size}

\item{Semi}{A list containing the semi-supervised estimates and their influence values}

\item{KM}{A list containing the Kaplan-Meier estimates and their influence values}
}
\value{
A numeric vector of SCANER estimates for survival probabilities at each time point
}
\description{
This function computes the combined SCANER estimator for survival probabilities
at each specified time point using both the semi-supervised and supervised estimators.
}
