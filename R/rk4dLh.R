#' Log-Likelihood Contributions for the rK4D Model
#'
#' Computes the observation-wise log-likelihood contributions for the
#' r-largest four-parameter kappa distribution (rK4D) model.
#'
#' @param data A numeric vector, matrix, or data frame of observations.
#'   If a vector is supplied, it is treated as a one-column matrix.
#'   If a matrix or data frame is supplied, each row is treated as one
#'   observation and columns represent decreasing order statistics.
#' @param par A numeric vector of length 4 giving the location, scale,
#'   first shape, and second shape parameters.
#'
#' @return A numeric vector of log-likelihood contributions for each row
#'   of \code{data}. If invalid parameter combinations occur, the function
#'   returns a large penalty value.
#'
#' @export
rk4dLh <- function(data, par) {

  tol <- .Machine$double.eps^0.5

  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  } else {
    data <- as.matrix(data)
  }

  if (!is.numeric(data)) {
    return(rep(1e6, nrow(data)))
  }

  if (!is.numeric(par) || length(par) != 4 || any(!is.finite(par))) {
    return(rep(1e6, nrow(data)))
  }

  R  <- ncol(data)
  nr <- nrow(data)

  mu <- par[1]
  sc <- par[2]
  xi <- par[3]
  h  <- par[4]

  if (!is.finite(sc) || sc <= 0) {
    return(rep(1e6, nr))
  }

  if (!is.finite(xi) || !is.finite(h) || abs(xi) < tol || abs(h) < tol) {
    return(rep(1e6, nr))
  }

  ri <- R - seq_len(R)
  cr <- 1 - ri * h

  if (any(cr <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  y <- 1 - xi * ((data - mu) / sc)

  if (any(y <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  yr <- 1 - xi * ((data[, R] - mu) / sc)

  if (any(yr <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  f <- 1 - h * yr^(1 / xi)

  if (any(!is.finite(f), na.rm = TRUE) || any(f <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  if (any(f^(1 / h) > 1, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  y <- log(sc) + (1 - 1 / xi) * log(y) - log(cr)
  y <- rowSums(y, na.rm = TRUE)

  log.den <- ((R * h - 1) / h) * log(f) + y

  if (any(!is.finite(log.den), na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  -log.den
}
