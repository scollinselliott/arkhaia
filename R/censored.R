#' Least Squares Fit of Poisson Distribution for Random Right-Censored Data
#'
#' For a matrix of cross-tabulated counts of observations which constitute a minimum threshold, this function esimates the rate parameter column-wise, either retaining or omitting zeros, by a least-squares approach.
#' 
#' @param x A matrix of cross-tabulated counts.
#' @param lambda_grid The resolutation at which to sample for the rate parameter. Default is \code{seq(0.1, 100, by = 0.01)}.
#' @param omit_zero Whether to omit zeros. Default is \code{TRUE}.
#' @examples 
#' x1 <- c(1,2,2,5,7,0,0)
#' x2 <- c(9,2,5,15,7,90,0)
#' x <- matrix(c(x1,x2), ncol = 2)
#' 
#' pois_rcens(x)
#' pois_rcens(x, omit_zero = FALSE)
#'
#' @returns A vector of the rate parameters for each column.
#' 
#' @importFrom Rdpack reprompt
#' @export
pois_rcens <- function(x, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE) {
    UseMethod("pois_rcens")
}

#' @rdname pois_rcens
#' @export
pois_rcens.matrix <- function(x, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE) {
  if (omit_zero == TRUE) {
    omit <- 1
  } else {
    omit <- 0
  }
  out <- lambda_col_arma(x, lambda_grid, omit)
  out <- as.vector(out)
  names(out) <- colnames(x)
  return(out)
}

#' @rdname pois_rcens
#' @export
pois_rcens.data.frame <- function(x, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE) {
  x <- as.matrix(x)
  out <- pois_rcens.matrix(x, lambda_grid = lambda_grid, omit_zero = omit_zero)
  return(out)
}

#' @rdname pois_rcens
#' @export
pois_rcens.xtabs <- function(x, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE) {
  x <- as.matrix(x)
  out <- pois_rcens.matrix(x, lambda_grid = lambda_grid, omit_zero = omit_zero)
  return(out)
}

#' @rdname pois_rcens
#' @export
pois_rcens.table <- function(x, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE) {
  x <- as.matrix(x)
  out <- pois_rcens.matrix(x, lambda_grid = lambda_grid, omit_zero = omit_zero)
  return(out)
}




#' Resampled Contingency Table via a Truncated Poisson for Random Right-Censored Data 
#'
#' For a matrix of cross-tabulated counts of observations which constitute a minimum threshold, returns a contingency table whose counts are sampled according to a truncated Poisson distribution, whose rate parameter is determined column-wise (see \code{\link[arkhaia]{pois_rcens}}). 
#' 
#' @param x A matrix of cross-tabulated counts.
#' @param lambda_grid The resolutation at which to sample for the rate parameter. Default is \code{seq(0.1, 100, by = 0.01)}.
#' @param omit_zero Whether to omit zeros. Default is \code{TRUE}.
#' 
#' @examples 
#' x1 <- c(1,2,2,5,7,0,0)
#' x2 <- c(9,2,5,15,7,90,0)
#' x <- matrix(c(x1,x2), ncol = 2)
#' 
#' pois_rcens(x)
#' pois_rcens(x, omit_zero = FALSE)
#'
#' @returns A contingency table \eqn{Y} of the same size as \eqn{X}, with \eqn{y_{ij}} drawn according to a truncated Poisson distribution that ensures \eqn{y_{ij} \geq x_{ij}}.
#' 
#' @importFrom Rdpack reprompt
#' @export
trunc_pois <- function(x, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE) {
    UseMethod("trunc_pois")
}

#' @rdname trunc_pois
#' @export
trunc_pois.matrix <- function(x, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE) {
  if (omit_zero == TRUE) {
    omit <- 1
  } else {
    omit <- 0
  }
  out <-trunc_pois_mat_arma(x, lambda_grid, omit)
  return(out)
}

#' @rdname trunc_pois
#' @export
trunc_pois.data.frame <- function(x, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE) {
  x <- as.matrix(x)
  out <- trunc_pois.matrix(x, lambda_grid = lambda_grid, omit_zero = omit_zero)
  return(out)
}

#' @rdname trunc_pois
#' @export
trunc_pois.xtabs <- function(x, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE) {
  x <- as.matrix(x)
  out <- trunc_pois.matrix(x, lambda_grid = lambda_grid, omit_zero = omit_zero)
  return(out)
}

#' @rdname trunc_pois
#' @export
trunc_pois.table <- function(x, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE) {
  x <- as.matrix(x)
  out <- trunc_pois.matrix(x, lambda_grid = lambda_grid, omit_zero = omit_zero)
  return(out)
}






#' Bias-Corrected Cramer's V for Random Right-Censored Data 
#'
#' Given a matrix of cross-tabulated counts of observations which constitute a minimum threshold, returns the distribution of bias-corrected Cramer's \eqn{V} by resampling from a truncated Poisson distribution, with rates determined column-wise (see \code{\link[arkhaia]{VB}}, \code{\link[arkhaia]{trunc_pois}}). 
#' 
#' @param x A matrix of cross-tabulated counts.
#' @param lambda The value of lambda (default is 2/3) for the Cressie-Read power divergence statistic used to estimate the bias-corrected Cramer's \eqn{V}. Default is \code{2/3}.
#' @param lambda_grid The resolutation at which to sample for the rate parameter. Default is \code{seq(0.1, 100, by = 0.01)}.
#' @param omit_zero Whether to omit zeros. Default is \code{TRUE}.
#' @param n_iter Number of samples of \eqn{V} to take. Default is \code{10^5}.
#' @examples 
#' x1 <- c(1,2,2,5,7,0,0)
#' x2 <- c(9,2,5,15,7,90,0)
#' x <- matrix(c(x1,x2), ncol = 2)
#' 
#' VB_trunc_pois(x, n_iter = 10^2)
#' VB_trunc_pois(x, omit_zero = FALSE, n_iter = 10^2)
#'
#' @returns A contingency table \eqn{Y} of the same size as \eqn{X}, with \eqn{y_{ij}} drawn according to a truncated Poisson distribution, \eqn{y_{ij} \geq x_{ij}}.
#' 
#' @importFrom Rdpack reprompt
#' @export
VB_trunc_pois <- function(x, lambda = 2/3, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE, n_iter = 10^5)  {
    UseMethod("VB_trunc_pois")
}

#' @rdname VB_trunc_pois
#' @export
VB_trunc_pois.matrix <- function(x, lambda = 2/3, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE, n_iter = 10^5) {
  if (omit_zero == TRUE) {
    omit <- 1
  } else {
    omit <- 0
  }
  out <- MC_pois_arma(x, lambda, lambda_grid, omit, n_iter)
  return(out)
}

#' @rdname VB_trunc_pois
#' @export
VB_trunc_pois.data.frame <- function(x, lambda = 2/3, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE, n_iter = 10^5)  {
  x <- as.matrix(x)
  out <- VB_trunc_pois.matrix(x, lambda, lambda_grid, omit_zero, n_iter)
  return(out)
}

#' @rdname VB_trunc_pois
#' @export
VB_trunc_pois.xtabs <- function(x, lambda = 2/3, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE, n_iter = 10^5) {
  x <- as.matrix(x)
  out <- VB_trunc_pois.matrix(x, lambda, lambda_grid, omit_zero, n_iter)
  return(out)
}

#' @rdname VB_trunc_pois
#' @export
VB_trunc_pois.table <- function(x, lambda = 2/3, lambda_grid = seq(0.01, 100, by = 0.01), omit_zero = TRUE, n_iter = 10^5)  {
  x <- as.matrix(x)
  out <- VB_trunc_pois.matrix(x, lambda, lambda_grid, omit_zero, n_iter)
  return(out)
}





