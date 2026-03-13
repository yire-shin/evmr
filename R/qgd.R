#' Quantile Function of the Gumbel Distribution
#'
#' Computes the quantiles of the Gumbel distribution with location
#' parameter \code{loc} and scale parameter \code{scale}.
#'
#' @param p A numeric vector of probabilities in \eqn{(0,1)}.
#' @param loc A numeric value specifying the location parameter.
#' @param scale A positive numeric value specifying the scale parameter.
#'
#' @return A numeric vector of quantiles corresponding to \code{p}.
#'
#' @details
#' The quantile function of the Gumbel distribution is
#' \deqn{Q(p) = \mu - \sigma \log(-\log(p)),}
#' where \eqn{\mu} is the location parameter and \eqn{\sigma > 0}
#' is the scale parameter.
#'
#' @export
#'
#' @examples
#' qgd(0.5, loc = 0, scale = 1)
#' qgd(c(0.1, 0.5, 0.9), loc = 0, scale = 1)
qgd <- function(p, loc = 0, scale = 1) {
  if (!is.numeric(p)) {
    stop("'p' must be numeric.")
  }

  if (any(!is.finite(p))) {
    stop("'p' must contain only finite values.")
  }

  if (any(p <= 0 | p >= 1)) {
    stop("'p' must be in (0,1).")
  }

  if (!is.numeric(loc) || any(!is.finite(loc))) {
    stop("'loc' must be numeric and finite.")
  }

  if (!is.numeric(scale) || any(!is.finite(scale))) {
    stop("'scale' must be numeric and finite.")
  }

  if (any(scale <= 0)) {
    stop("'scale' must be positive.")
  }

  mu <- loc
  sig <- scale

  yp <- -log(p)
  z <- mu - sig * log(yp)

  return(z)
}
