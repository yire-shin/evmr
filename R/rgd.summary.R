#' @name rgd.summary
#' @aliases rgd.summary
#' @title Summary of rgd.fit
#' @description
#' Summarize the fitted Gumbel distribution models for r-largest order
#' statistics over \eqn{r = 1, \dots, R}. For each value of \code{r},
#' the function fits the model using \code{rgd.fit()} and computes
#' return levels using \code{rgd.rl()}.
#'
#' @param data A numeric vector, matrix, or data frame containing the
#' r-largest order statistics. Each row should contain decreasing order
#' statistics for one block or time period.
#' @param r Optional integer specifying the maximum number of order
#' statistics to summarize. If \code{NULL}, all available columns are used.
#' @param ydat A matrix of covariates for generalized linear modelling of the
#' parameters, or \code{NULL} for stationary fitting.
#' @param mul,sigl Numeric vectors of integers giving the columns of
#' \code{ydat} to be used for the location and scale parameters, respectively.
#' @param mulink,siglink Inverse link functions for the location and scale
#' parameters, respectively.
#' @param num_inits Number of initial parameter sets used in optimization.
#' @param muinit,siginit Optional initial values for the location and scale
#' parameters.
#' @param show Logical; if \code{TRUE}, print details from model fitting.
#' @param method Optimization method passed to \code{\link{optim}}.
#' @param maxit Maximum number of iterations for optimization.
#' @param \dots Additional arguments passed to \code{rgd.fit()}.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{r}: number of order statistics used
#'   \item \code{nllh}: negative log-likelihood
#'   \item \code{mu}, \code{sigma}: parameter estimates
#'   \item \code{mu.se}, \code{sigma.se}: standard errors
#'   \item return levels and their standard errors for 20-, 50-, 100-, and 200-year periods
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(bangkok)
#' rgd.summary(bangkok)
#' }
rgd.summary <- function(data, r= NULL, ydat = NULL, mul = NULL, sigl = NULL,
                        mulink = identity, siglink = identity, num_inits = 100,
                        muinit = NULL, siginit = NULL, show = FALSE,
                        method = "Nelder-Mead", maxit = 10000, ...) {

  R <- ncol(data)

  result <- matrix(NA_real_, R, 14)

  for (i in seq_len(R)) {

    xsub <- data[, 1:i, drop = FALSE]

    fit <- rgd.fit(
      xdat = xsub,
      r = r,
      ydat = ydat,
      mul = mul,
      sigl = sigl,
      mulink = mulink,
      siglink = siglink,
      num_inits = num_inits,
      muinit = muinit,
      siginit = siginit,
      show = show,
      method = method,
      maxit = maxit,
      ...
    )

    rl <- rgd.rl(fit, show = FALSE)

    result[i, 1]  <- i
    result[i, 2]  <- fit$nllh
    result[i, 3:4] <- fit$mle
    result[i, 5:6] <- fit$se
    result[i, 7]  <- rl$rl[1]
    result[i, 8]  <- rl$rlse[1]
    result[i, 9]  <- rl$rl[2]
    result[i,10]  <- rl$rlse[2]
    result[i,11]  <- rl$rl[3]
    result[i,12]  <- rl$rlse[3]
    result[i,13]  <- rl$rl[4]
    result[i,14]  <- rl$rlse[4]
  }

  colnames(result) <- c(
    "r", "nllh", "mu", "sigma", "mu.se", "sigma.se",
    "rl20", "rl20.se", "rl50", "rl50.se",
    "rl100", "rl100.se", "rl200", "rl200.se"
  )

  result <- as.data.frame(result)
  result[, 2:ncol(result)] <- round(result[, 2:ncol(result)], 3)

  result
}
