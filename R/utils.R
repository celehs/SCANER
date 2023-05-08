#' Generate random event times
#'
#' This function generates random event times using an exponential distribution.
#'
#' @param nn Integer specifying the number of random event times to generate
#' @return A vector of nn random event times
#' @export
rT.fun = function(nn) {
  rexp(nn)
}

#' Generate random censoring times
#'
#' This function generates random censoring times using an exponential distribution.
#'
#' @param nn Integer specifying the number of random censoring times to generate
#' @return A vector of nn random censoring times
#' @export
rC.fun = function(nn) {
  rexp(nn)
}

#' Create a matrix with identical rows
#'
#' This function creates a matrix with dm rows and n columns, where n is the length of the input vector vc. Each row of the matrix is identical to vc.
#'
#' @param vc Numeric vector to be used as the template for each row of the matrix
#' @param dm Integer specifying the number of rows in the output matrix
#' @return An dm x n matrix with each row equal to vc
#' @export
#'
#' @examples
#' # Example usage of VTM function
#' vec <- c(1, 2, 3)
#' mat <- VTM(vec, 4)
VTM <- function(vc, dm) {
  matrix(vc, ncol = length(vc), nrow = dm, byrow = T)
}
