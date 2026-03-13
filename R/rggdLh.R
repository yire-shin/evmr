#' Negative Log-Likelihood for the rGGD Model
#'
#' Computes the negative log-likelihood for the r-largest generalized
#' Gumbel distribution (rGGD) model.
#'
#' @param data A numeric vector, matrix, or data frame of observations.
#'   If a vector is supplied, it is treated as a one-column matrix.
#'   If a matrix or data frame is supplied, each row is treated as one
#'   observation and columns represent decreasing order statistics.
#' @param par A numeric vector of length 3 giving the location, scale,
#'   and shape parameters, respectively.
#'
#' @return A single numeric value giving the negative log-likelihood.
#'   If the parameter combination is invalid, the function returns \code{Inf}.
#'
#' @details
#' This function is intended for internal likelihood evaluation in optimization.
#' Invalid parameter combinations return \code{Inf} rather than stopping with
#' an error, which makes the function more robust when used inside optimizers
#' such as \code{\link{optim}}.
#'
#'#' @references
#' Shin, Y., & Park, J.-S. (2025).
#' Generalized Gumbel model for r-largest order statistics with application
#' to peak streamflow.
#' \emph{Scientific Reports}.
#' \doi{10.1038/s41598-024-83273-y}
#'
#' @export
rggdLh <- function(data, par) {

  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
  } else {
    data <- as.matrix(data)
  }

  nr <- nrow(data)
  R <- ncol(data)

  if (!is.numeric(data)) {
    return(rep(1e6, nr))
  }

  if (!is.numeric(par) || length(par) != 3 || any(!is.finite(par))) {
    return(rep(1e6, nr))
  }

  mu <- par[1]
  sc <- par[2]
  h  <- par[3]

  # invalid scale
  if (!is.finite(sc) || sc <= 0) {
    return(rep(1e6, nr))
  }

  # avoid numerical instability near h = 0
  if (!is.finite(h) || abs(h) < 1e-8) {
    return(rep(1e6, nr))
  }

  # constraint: 1 - (r - i)h > 0
  ri <- R - seq_len(R)
  cr <- 1 - ri * h

  if (any(cr <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  y <- exp(-(data - mu) / sc)
  f <- 1 - h * exp(-(data[, R] - mu) / sc)

  if (any(!is.finite(y), na.rm = TRUE) ||
      any(!is.finite(f), na.rm = TRUE) ||
      any(y <= 0, na.rm = TRUE) ||
      any(f <= 0, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  if (any(f^(1 / h) > 1, na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  y <- log(sc) - log(y) - log(cr)
  y <- rowSums(y, na.rm = TRUE)

  log.den3 <- ((R * h - 1) / h) * log(f) + y

  if (any(!is.finite(log.den3), na.rm = TRUE)) {
    return(rep(1e6, nr))
  }

  -log.den3
}
