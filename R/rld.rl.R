#' Return Levels for the Logistic Distribution
#'
#' Computes return levels and their standard errors for a stationary
#' logistic model fitted by \code{\link{rld.fit}}.
#'
#' @param z An object returned by \code{\link{rld.fit}}. The fitted model
#'   should represent a stationary model.
#' @param year A numeric vector of return periods for which return levels
#'   are to be computed.
#' @param show Logical. If \code{TRUE}, the estimated return levels and
#'   their standard errors are printed.
#'
#' @return The input object \code{z} with two additional components:
#' \itemize{
#'   \item \code{rl}: a numeric vector of estimated return levels
#'   \item \code{rlse}: a numeric vector of standard errors of the estimated return levels
#' }
#'
#' @details
#' For a return period \eqn{T}, the return level is defined as the quantile
#' exceeded with probability \eqn{1/T}. Under the logistic distribution,
#' the return level is
#' \deqn{x_T = \mu + \sigma \log\left(\frac{1}{\exp(-\log(1-1/T)) - 1}\right),}
#' and standard errors are obtained using the delta method.
#'
#' @seealso \code{\link{rld.fit}}, \code{\link{rld.prof}}
#' @export
#'
#' @examples
#' x <- rldr(n = 50, r = 2, loc = 10, scale = 2)
#' fit <- rld.fit(x$rmat, num_inits = 5)
#' out <- rld.rl(fit, year = c(20, 50, 100))
#' out$rl
#' out$rlse
rld.rl <- function(z, year = c(20, 50, 100, 200), show = FALSE) {

  if (!inherits(z, "rld.fit")) {
    warning("'z' does not inherit from class 'rld.fit'.")
  }

  if (isTRUE(z$trans)) {
    stop("Return levels are only implemented for stationary models.")
  }

  if (!is.numeric(year) || any(!is.finite(year)) || any(year <= 1)) {
    stop("'year' must be a numeric vector with values greater than 1.")
  }

  if (is.null(z$mle) || length(z$mle) < 2) {
    stop("'z$mle' must contain location and scale estimates.")
  }

  if (is.null(z$cov) || all(is.na(z$cov))) {
    stop("Covariance matrix is not available in 'z$cov'.")
  }

  mu <- z$mle[1]
  sig <- z$mle[2]

  if (!is.finite(sig) || sig <= 0) {
    stop("The fitted scale parameter must be positive.")
  }

  del <- matrix(nrow = 2, ncol = length(year))

  f <- 1 - (1 / year)

  if (any(!is.finite(f)) || any(f <= 0) || any(f >= 1)) {
    stop("Invalid return period values in 'year'.")
  }

  del[1, ] <- 1
  del[2, ] <- log(1 / expm1(-log(f)))

  del.t <- t(del)

  z$rl <- mu + sig * log(1 / expm1(-log(f)))
  z$rlse <- sqrt(diag(del.t %*% z$cov %*% del))

  names(z$rl) <- paste0(as.character(year), "y")
  names(z$rlse) <- paste0(as.character(year), "y")

  if (show) {
    print(z[c("rl", "rlse")])
  }

  invisible(z)
}
