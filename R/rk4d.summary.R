#' Summary of Fitted rK4D Models over Different Values of r
#'
#' Summarizes fitted four-parameter kappa distribution models for
#' r-largest order statistics over \eqn{r = 1, \dots, R}. For each value
#' of \code{r}, the function fits the model using \code{\link{rk4d.fit}}
#' and computes return levels using \code{\link{rk4d.rl}}.
#'
#' @param data A numeric vector, matrix, or data frame containing the
#'   r-largest order statistics. Each row should contain decreasing order
#'   statistics for one block or time period.
#' @param r Optional integer giving the maximum number of order statistics
#'   to summarize. If \code{NULL}, all available columns are used.
#' @param penk Penalty function for the \code{xi} parameter in maximum
#'   penalized likelihood estimation.
#' @param penh Penalty function for the \code{h} parameter in maximum
#'   penalized likelihood estimation.
#' @param ydat A matrix or data frame of covariates for generalized linear
#'   modelling of the parameters, or \code{NULL} for stationary fitting.
#' @param mul,sigl,shl,hl Integer vectors indicating which columns of
#'   \code{ydat} are used for the location, scale, first shape, and second
#'   shape parameters, respectively.
#' @param mulink,siglink,shlink,hlink Inverse link functions for the location,
#'   scale, first shape, and second shape parameters, respectively.
#' @param num_inits Number of initial parameter sets used in optimization.
#' @param muinit,siginit,shinit,hinit Optional initial values for the location,
#'   scale, first shape, and second shape parameters.
#' @param show Logical. If \code{TRUE}, print details from model fitting.
#' @param method Optimization method passed to \code{\link{optim}}.
#' @param maxit Maximum number of iterations for optimization.
#' @param ... Additional arguments passed to \code{\link{rk4d.fit}}.
#'
#' @return A data frame containing:
#' \itemize{
#'   \item \code{r}: number of order statistics used
#'   \item \code{nllh}: negative log-likelihood
#'   \item \code{mu}, \code{sigma}, \code{xi}, \code{h}: parameter estimates
#'   \item \code{mu.se}, \code{sigma.se}, \code{xi.se}, \code{h.se}: standard errors
#'   \item \code{rl20}, \code{rl50}, \code{rl100}, \code{rl200}: return levels
#'   \item \code{rl20.se}, \code{rl50.se}, \code{rl100.se}, \code{rl200.se}:
#'         standard errors of return levels
#' }
#'
#' @export
#'
#' @examples
#' x <- rk4dr(n = 50, r = 3, loc = 10, scale = 2, shape1 = 0.1, shape2 = 0.1)
#' rk4d.summary(x$rmat, num_inits = 5)
#' rk4d.summary(x$rmat, penk = "CD", penh = "MS", num_inits = 5)
rk4d.summary <- function(data, r = NULL, penk = NULL, penh = NULL,
                         ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, hl = NULL,
                         mulink = identity, siglink = identity,
                         shlink = identity, hlink = identity,
                         num_inits = 100, muinit = NULL, siginit = NULL,
                         shinit = NULL, hinit = NULL, show = FALSE,
                         method = "Nelder-Mead", maxit = 10000, ...) {

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
      if (!is.numeric(r) || length(r) != 1 || r < 1 || r != as.integer(r)) {
        stop("'r' must be a positive integer.")
      }
      R <- min(as.integer(r), ncol(data))
    }
  }

  result <- matrix(NA_real_, nrow = R, ncol = 18)

  for (i in seq_len(R)) {
    xsub <- data[, 1:i, drop = FALSE]

    fit <- tryCatch(
      rk4d.fit(
        xdat = xsub,
        r = i,
        penk = penk,
        penh = penh,
        ydat = ydat,
        mul = mul,
        sigl = sigl,
        shl = shl,
        hl = hl,
        mulink = mulink,
        siglink = siglink,
        shlink = shlink,
        hlink = hlink,
        num_inits = num_inits,
        muinit = muinit,
        siginit = siginit,
        shinit = shinit,
        hinit = hinit,
        show = show,
        method = method,
        maxit = maxit,
        ...
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) {
      next
    }

    rl <- tryCatch(
      rk4d.rl(fit, show = FALSE),
      error = function(e) NULL
    )

    result[i, 1] <- i
    result[i, 2] <- fit$nllh

    if (!is.null(fit$mle) && length(fit$mle) >= 4) {
      result[i, 3:6] <- fit$mle[1:4]
    }

    if (!all(is.na(fit$se)) && length(fit$se) >= 4) {
      result[i, 7:10] <- fit$se[1:4]
    }

    if (!is.null(rl)) {
      result[i, 11] <- rl$rl[1]
      result[i, 12] <- rl$rlse[1]
      result[i, 13] <- rl$rl[2]
      result[i, 14] <- rl$rlse[2]
      result[i, 15] <- rl$rl[3]
      result[i, 16] <- rl$rlse[3]
      result[i, 17] <- rl$rl[4]
      result[i, 18] <- rl$rlse[4]
    }
  }

  colnames(result) <- c(
    "r", "nllh", "mu", "sigma", "xi", "h",
    "mu.se", "sigma.se", "xi.se", "h.se",
    "rl20", "rl20.se", "rl50", "rl50.se",
    "rl100", "rl100.se", "rl200", "rl200.se"
  )

  result <- as.data.frame(result)
  result[, -1] <- round(result[, -1], 3)

  result
}
