#' Quantile Function of the Generalized Gumbel Distribution
#'
#' Computes the quantiles of the generalized Gumbel distribution
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
#' \deqn{Q(p) = \mu - \sigma \log \left( \frac{1 - p^h}{h} \right), \quad h \neq 0,}
#' with the limiting case
#' \deqn{Q(p) = \mu - \sigma \log(-\log p), \quad h = 0,}
#' where \eqn{\mu} is the location parameter, \eqn{\sigma > 0} is the
#' scale parameter, and \eqn{h} is the shape parameter.
#'
#' @references
#' Jeong, B.-Y., Murshed, M. S., Seo, Y. A., and Park, J.-S. (2014).
#' A three-parameter kappa distribution with hydrologic application:
#' a generalized Gumbel distribution.
#' \emph{Stochastic Environmental Research and Risk Assessment},
#' 28(8), 2063--2074.
#'
#' @export
#'
#' @examples
#' qggd(0.5, loc = 0, scale = 1, shape = 0.1)
#' qggd(c(0.1, 0.5, 0.9), loc = 0, scale = 1, shape = 0.1)
qggd <- function(p, loc = 0, scale = 1, shape = 0) {
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

  # h = 0: Gumbel limit
  idx0 <- abs(shape) < .Machine$double.eps^0.5
  if (any(idx0)) {
    z[idx0] <- loc[idx0] - scale[idx0] * log(-log(p[idx0]))
  }

  # h != 0: generalized case
  idx1 <- !idx0
  if (any(idx1)) {
    z[idx1] <- loc[idx1] - scale[idx1] *
      log((1 - p[idx1]^shape[idx1]) / shape[idx1])
  }

  z
}
