#' Log-Likelihood Contributions for the rGLO Model
#'
#' Computes the observation-wise log-likelihood contributions for the
#' r-largest generalized logistic distribution (rGLO) model.
#'
#' @param data A numeric vector, matrix, or data frame of observations.
#'   If a vector is supplied, it is treated as a one-column matrix.
#'   If a matrix or data frame is supplied, each row is treated as one
#'   observation and columns represent decreasing order statistics.
#' @param par A numeric vector of length 3 giving the location, scale,
#'   and shape parameters, respectively.
#'
#' @return A numeric vector of log-likelihood contributions, one for each row
#'   of \code{data}. If the parameter combination is invalid, the function
#'   returns \code{Inf}.
#'
#' @details
#' This function is mainly intended for internal likelihood evaluation.
#' Invalid parameter combinations return \code{Inf}, which is often more
#' robust than stopping with an error when used inside iterative procedures.
#'
#' @export
rgloLh <- function(data, par) {

  tol <- .Machine$double.eps^0.5

  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  } else {
    data <- as.matrix(data)
  }

  nr <- nrow(data)

  if (!is.numeric(data)) {
    return(rep(1e6, nr))
  }

  if (!is.numeric(par) || length(par) != 3 || any(!is.finite(par))) {
    return(rep(1e6, nr))
  }

  R <- ncol(data)

  mu <- par[1]
  sc <- par[2]
  xi <- par[3]

  if (!is.finite(sc) || sc <= 0) {
    return(rep(1e6, nr))
  }

  if (!is.finite(xi) || abs(xi) < tol) {
    return(rep(1e6, nr))
  }

  ri <- R - seq_len(R)
  cr <- 1 + ri

  if (any(cr <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  y <- 1 - xi * ((data - mu) / sc)
  yr <- 1 - xi * ((data[, R] - mu) / sc)

  if (any(y <= 0, na.rm = TRUE) || any(yr <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  f <- 1 + yr^(1 / xi)

  if (any(!is.finite(f), na.rm = TRUE) || any(f <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  log.den <- -R * log(sc) +
    sum(log(cr)) +
    (1 + R) * log(1 / (1 + yr^(1 / xi))) +
    rowSums(((1 / xi) - 1) * log(y), na.rm = TRUE)

  if (any(!is.finite(log.den), na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  log.den
}
