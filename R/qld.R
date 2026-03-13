#' Quantile Function of the Logistic Distribution
#'
#' Computes the quantiles of the logistic distribution with location
#' parameter \code{loc} and scale parameter \code{scale}.
#'
#' @param p A numeric vector of probabilities in \eqn{(0,1)}.
#' @param loc A numeric value specifying the location parameter.
#' @param scale A positive numeric value specifying the scale parameter.
#'
#' @return A numeric vector of quantiles corresponding to \code{p}.
#'
#' @details
#' The quantile function of the logistic distribution is
#' \deqn{Q(p) = \mu + \sigma \log\left(\frac{p}{1-p}\right),}
#' where \eqn{\mu} is the location parameter and \eqn{\sigma > 0}
#' is the scale parameter.
#'
#' @export
#'
#' @examples
#' qld(0.5, loc = 0, scale = 1)
#' qld(c(0.1, 0.5, 0.9), loc = 0, scale = 1)
qld <- function(p, loc = 0, scale = 1) {
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

  lens <- c(length(p), length(loc), length(scale))
  L <- max(lens)

  ok_len <- function(x) length(x) %in% c(1, L)
  if (!ok_len(p) || !ok_len(loc) || !ok_len(scale)) {
    stop("'p', 'loc', and 'scale' must have length 1 or a common length.")
  }

  p <- rep_len(p, L)
  loc <- rep_len(loc, L)
  scale <- rep_len(scale, L)

  z <- loc + scale * log(p / (1 - p))

  z
}
