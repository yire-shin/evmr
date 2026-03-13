#' Summary of Fitted rGGD Models over Different Values of r
#'
#' Summarizes fitted generalized Gumbel distribution models for
#' r-largest order statistics over \eqn{r = 1, \dots, R}. For each value
#' of \code{r}, the function fits the model using \code{\link{rggd.fit}}
#' and computes return levels using \code{\link{rggd.rl}}.
#'
#' @param data A numeric vector, matrix, or data frame containing the
#'   r-largest order statistics. Each row should contain decreasing order
#'   statistics for one block or time period.
#' @param r Optional integer giving the maximum number of order statistics
#'   to summarize. If \code{NULL}, all available columns are used.
#' @param ydat A matrix or data frame of covariates for generalized linear
#'   modelling of the parameters, or \code{NULL} for stationary fitting.
#' @param mul,sigl,hl Integer vectors indicating which columns of
#'   \code{ydat} are used for the location, scale, and shape parameters,
#'   respectively.
#' @param mulink,siglink,hlink Inverse link functions for the location,
#'   scale, and shape parameters, respectively.
#' @param num_inits Number of initial parameter sets used in optimization.
#' @param muinit,siginit,hinit Optional initial values for the location,
#'   scale, and shape parameters.
#' @param show Logical. If \code{TRUE}, print details from model fitting.
#' @param method Optimization method passed to \code{\link{optim}}.
#' @param maxit Maximum number of iterations for optimization.
#' @param ... Additional arguments passed to \code{\link{rggd.fit}}.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{r}: number of order statistics used
#'   \item \code{nllh}: negative log-likelihood
#'   \item \code{mu}, \code{sigma}, \code{h}: parameter estimates
#'   \item \code{mu.se}, \code{sigma.se}, \code{h.se}: standard errors
#'   \item \code{rl20}, \code{rl50}, \code{rl100}, \code{rl200}: return levels
#'   \item \code{rl20.se}, \code{rl50.se}, \code{rl100.se}, \code{rl200.se}:
#'         standard errors of return levels
#' }
#'
#' @export
#'
#' @examples
#' x <- rggdr(n = 50, r = 3, loc = 10, scale = 2, shape = 0.1)
#' rggd.summary(x$rmat)
rggd.summary <- function(data, r = NULL, ydat = NULL, mul = NULL, sigl = NULL, hl = NULL,
                         mulink = identity, siglink = identity, hlink = identity,
                         num_inits = 100, muinit = NULL, siginit = NULL, hinit = NULL,
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

  result <- matrix(NA_real_, nrow = R, ncol = 16)

  for (i in seq_len(R)) {
    xsub <- data[, 1:i, drop = FALSE]

    fit <- rggd.fit(
      xdat = xsub,
      r = i,
      ydat = ydat,
      mul = mul,
      sigl = sigl,
      hl = hl,
      mulink = mulink,
      siglink = siglink,
      hlink = hlink,
      num_inits = num_inits,
      muinit = muinit,
      siginit = siginit,
      hinit = hinit,
      show = show,
      method = method,
      maxit = maxit,
      ...
    )

    rl <- rggd.rl(fit, show = FALSE)

    result[i, 1] <- i
    result[i, 2] <- fit$nllh
    result[i, 3:5] <- fit$mle
    result[i, 6:8] <- fit$se
    result[i, 9] <- rl$rl[1]
    result[i, 10] <- rl$rlse[1]
    result[i, 11] <- rl$rl[2]
    result[i, 12] <- rl$rlse[2]
    result[i, 13] <- rl$rl[3]
    result[i, 14] <- rl$rlse[3]
    result[i, 15] <- rl$rl[4]
    result[i, 16] <- rl$rlse[4]
  }

  colnames(result) <- c(
    "r", "nllh", "mu", "sigma", "h", "mu.se", "sigma.se", "h.se",
    "rl20", "rl20.se", "rl50", "rl50.se",
    "rl100", "rl100.se", "rl200", "rl200.se"
  )

  result <- as.data.frame(result)
  result[, -1] <- round(result[, -1], 3)

  result
}
