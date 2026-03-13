#' Return Levels for the Gumbel Distribution
#'
#' Computes return levels and their standard errors for a stationary
#' Gumbel model fitted by \code{\link{rgd.fit}}.
#'
#' @param z An object returned by \code{\link{rgd.fit}}. The fitted model
#'   should represent a stationary model.
#' @param year A numeric vector of return periods for which return levels
#'   are to be computed.
#' @param show Logical. If \code{TRUE}, the estimated return levels and
#'   their standard errors are printed.
#'
#' @return The input object \code{z} with two additional components:
#' \item{rl}{A numeric vector of estimated return levels.}
#' \item{rlse}{A numeric vector of standard errors of the estimated return levels.}
#'
#' @details
#' For a return period \eqn{T}, the return level is defined as the quantile
#' exceeded with probability \eqn{1/T}. Under the Gumbel distribution, the
#' return level is
#' \deqn{x_T = \mu - \sigma \log\{-\log(1 - 1/T)\}.}
#' Standard errors are obtained using the delta method.
#'
#' @seealso \code{\link{rgd.fit}}, \code{\link{rgd.prof}}
#' @export
#'
#' @examples
#' x <- rgdr(n = 50, r = 2, loc = 10, scale = 2)
#' fit <- rgd.fit(x$rmat)
#' out <- rgd.rl(fit, year = c(20, 50, 100))
#' out$rl
#' out$rlse
rgd.rl <- function(z, year = c(20, 50, 100, 200), show = FALSE) {

  if (!inherits(z, "rgd.fit")) {
    warning("'z' does not inherit from class 'rgd.fit'.")
  }

  if (isTRUE(z$trans)) {
    stop("Return levels are only implemented for stationary models.")
  }

  if (!is.numeric(year) || any(!is.finite(year)) || any(year <= 1)) {
    stop("'year' must be a numeric vector with values greater than 1.")
  }

  if (is.null(z$mle) || length(z$mle) < 2) {
    stop("'z$mle' must contain at least location and scale estimates.")
  }

  if (is.null(z$cov) || all(is.na(z$cov))) {
    stop("Covariance matrix is not available in 'z$cov'.")
  }

  del <- matrix(ncol = length(year), nrow = 2)

  mu <- z$mle[1]
  sig <- z$mle[2]

  p <- 1 / year
  yp <- -log(1 - p)

  del[1, ] <- 1
  del[2, ] <- -log(yp)

  del.t <- t(del)

  z$rl <- mu - sig * log(yp)
  z$rlse <- sqrt(diag(del.t %*% z$cov %*% del))

  names(z$rl) <- paste0(as.character(year), "y")
  names(z$rlse) <- paste0(as.character(year), "y")

  if (show) {
    print(z[c("rl", "rlse")])
  }

  invisible(z)
}
