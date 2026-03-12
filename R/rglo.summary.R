#' @name rglo.summary
#' @aliases rglo.summary
#' @title Summary of rglo.fit
#' @description
#' Summarize the fitted generalized logistic distribution models for
#' r-largest order statistics over \eqn{r = 1, \dots, R}. For each value
#' of \code{r}, the function fits the model using \code{rglo.fit()} and
#' computes return levels using \code{rglo.rl()}.
#'
#' @param data A numeric vector, matrix, or data frame containing the
#' r-largest order statistics. Each row should contain decreasing order
#' statistics for one block or time period.
#' @param r Optional integer specifying the maximum number of order
#' statistics to summarize. If \code{NULL}, all available columns are used.
#' @param ydat A matrix of covariates for generalized linear modelling of the
#' parameters, or \code{NULL} for stationary fitting.
#' @param mul,sigl,shl Numeric vectors of integers giving the columns of
#' \code{ydat} to be used for the location, scale, and shape parameters,
#' respectively.
#' @param mulink,siglink,shlink Inverse link functions for the location,
#' scale, and shape parameters, respectively.
#' @param num_inits Number of initial parameter sets used in optimization.
#' @param muinit,siginit,shinit Optional initial values for the location,
#' scale, and shape parameters.
#' @param show Logical; if \code{TRUE}, print details from model fitting.
#' @param method Optimization method passed to \code{\link{optim}}.
#' @param maxit Maximum number of iterations for optimization.
#' @param \dots Additional arguments passed to \code{rglo.fit()}.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{r}: number of order statistics used
#'   \item \code{nllh}: negative log-likelihood
#'   \item \code{mu}, \code{sigma}, \code{xi}: parameter estimates
#'   \item \code{mu.se}, \code{sigma.se}, \code{xi.se}: standard errors
#'   \item return levels and their standard errors for 20-, 50-, 100-, and 200-year periods
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rglo.summary(bangkok)
#' }
rglo.summary <- function(data, r = NULL, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL,
                         mulink = identity, siglink = identity, shlink = identity,
                         num_inits = 100, muinit = NULL, siginit = NULL, shinit = NULL,
                         show = FALSE, method = "Nelder-Mead", maxit = 10000, ...) {

  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }

  if (is.vector(data)) {
    data <- matrix(data, ncol = 1)
    R <- 1
  } else {
    data <- as.matrix(data)
    if (is.null(r)) {
      R <- ncol(data)
    } else {
      R <- min(r, ncol(data))
    }
  }

  result <- matrix(NA_real_, R, 16)

  for (i in seq_len(R)) {

    xsub <- as.matrix(data[, 1:i, drop = FALSE])

    fit <- rglo.fit(
      xdat = xsub,
      r = i,
      ydat = ydat,
      mul = mul,
      sigl = sigl,
      shl = shl,
      mulink = mulink,
      siglink = siglink,
      shlink = shlink,
      num_inits = num_inits,
      muinit = muinit,
      siginit = siginit,
      shinit = shinit,
      show = show,
      method = method,
      maxit = maxit,
      ...
    )

    rl <- rglo.rl(fit, show = FALSE)

    result[i, 1]  <- i
    result[i, 2]  <- fit$nllh
    result[i, 3:5] <- fit$mle
    result[i, 6:8] <- fit$se
    result[i, 9]  <- rl$rl[1]
    result[i,10]  <- rl$rlse[1]
    result[i,11]  <- rl$rl[2]
    result[i,12]  <- rl$rlse[2]
    result[i,13]  <- rl$rl[3]
    result[i,14]  <- rl$rlse[3]
    result[i,15]  <- rl$rl[4]
    result[i,16]  <- rl$rlse[4]
  }

  colnames(result) <- c(
    "r", "nllh", "mu", "sigma", "xi", "mu.se", "sigma.se", "xi.se",
    "rl20", "rl20.se", "rl50", "rl50.se",
    "rl100", "rl100.se", "rl200", "rl200.se"
  )

  result <- as.data.frame(result)
  result[, 2:ncol(result)] <- round(result[, 2:ncol(result)], 3)

  result
}
