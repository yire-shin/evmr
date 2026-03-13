#' Quantile Function of the Generalized Logistic Distribution
#'
#' Computes the quantiles of the generalized logistic distribution
#' with location parameter \code{loc}, scale parameter \code{scale},
#' and shape parameter \code{shape}.
#'
#' @param p A numeric vector of probabilities in \eqn{(0,1)}.
#' @param loc A numeric value specifying the location parameter.
#' @param scale A positive numeric value specifying the scale parameter.
#' @param shape A numeric value specifying the shape parameter.
#'
#' @return A numeric vector of quantiles corresponding to \code{p}.
#'
#' @details
#' The quantile function is computed as
#' \deqn{Q(p) = \mu + \frac{\sigma}{\xi}\left[1 - \left(\frac{1-p}{p}\right)^{\xi}\right], \quad \xi \neq 0,}
#' with the limiting case
#' \deqn{Q(p) = \mu - \sigma \log\left(\frac{1-p}{p}\right), \quad \xi = 0,}
#' where \eqn{\mu} is the location parameter, \eqn{\sigma > 0} is the
#' scale parameter, and \eqn{\xi} is the shape parameter.
#'
#' @references
#' Ahmad, M. I., Sinclair, C. D., and Werritty, A. (1988).
#' Log-logistic flood frequency analysis.
#' \emph{Journal of Hydrology}.
#' \doi{10.1016/0022-1694(88)90015-7}
#'
#' @export
#'
#' @examples
#' qglo(0.5, loc = 0, scale = 1, shape = 0.1)
#' qglo(c(0.1, 0.5, 0.9), loc = 0, scale = 1, shape = 0.1)
qglo <- function(p, loc = 0, scale = 1, shape = 0) {
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

  if (!is.numeric(shape) || any(!is.finite(shape))) {
    stop("'shape' must be numeric and finite.")
  }

  lens <- c(length(p), length(loc), length(scale), length(shape))
  L <- max(lens)

  ok_len <- function(x) length(x) %in% c(1, L)
  if (!ok_len(p) || !ok_len(loc) || !ok_len(scale) || !ok_len(shape)) {
    stop("'p', 'loc', 'scale', and 'shape' must have length 1 or a common length.")
  }

  p <- rep_len(p, L)
  loc <- rep_len(loc, L)
  scale <- rep_len(scale, L)
  shape <- rep_len(shape, L)

  z <- numeric(L)

  # xi = 0: logistic limit
  idx0 <- abs(shape) < .Machine$double.eps^0.5
  if (any(idx0)) {
    z[idx0] <- loc[idx0] - scale[idx0] * log((1 - p[idx0]) / p[idx0])
  }

  # xi != 0: generalized logistic case
  idx1 <- !idx0
  if (any(idx1)) {
    z[idx1] <- loc[idx1] +
      scale[idx1] / shape[idx1] -
      scale[idx1] / shape[idx1] *
      ((1 - p[idx1]) / p[idx1])^shape[idx1]
  }

  z
}
