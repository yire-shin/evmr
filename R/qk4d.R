#' Quantile Function of the Four-Parameter Kappa Distribution
#'
#' Computes the quantiles of the four-parameter kappa distribution
#' with location parameter \code{loc}, scale parameter \code{scale},
#' and shape parameters \code{shape1} and \code{shape2}.
#'
#' @param p A numeric vector of probabilities in \eqn{(0,1)}.
#' @param loc A numeric value specifying the location parameter.
#' @param scale A positive numeric value specifying the scale parameter.
#' @param shape1 A numeric value specifying the first shape parameter.
#' @param shape2 A numeric value specifying the second shape parameter.
#'
#' @return A numeric vector of quantiles corresponding to \code{p}.
#'
#' @details
#' The quantile function of the four-parameter kappa distribution is
#' \deqn{Q(p) = \mu + \frac{\sigma}{\xi}\left[1 - \left(\frac{1-p^h}{h}\right)^\xi \right],}
#' where \eqn{\mu} is the location parameter, \eqn{\sigma > 0} is the
#' scale parameter, and \eqn{\xi} and \eqn{h} are shape parameters.
#'
#' For numerical stability, the limiting cases \eqn{\xi = 0} and/or
#' \eqn{h = 0} are handled separately.
#'
#' @references
#' Shin, Y., and Park, J.-S.(2023).
#' Modeling climate extremes using the four-parameter kappa distribution
#' for r-largest order statistics.
#' \emph{Weather and Climate Extremes}.
#' \doi{10.1016/j.wace.2022.100533}
#'
#' Hosking, J. R. M. (1994).
#' \emph{The four-parameter Kappa distribution}.
#' Cambridge University Press.
#'
#' @export
#'
#' @examples
#' qk4d(0.5, loc = 0, scale = 1, shape1 = 0.1, shape2 = 0.1)
#' qk4d(c(0.1, 0.5, 0.9), loc = 0, scale = 1, shape1 = 0.1, shape2 = 0.1)
qk4d <- function(p, loc = 0, scale = 1, shape1 = 0.1, shape2 = 0.1) {
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

  if (!is.numeric(shape1) || any(!is.finite(shape1))) {
    stop("'shape1' must be numeric and finite.")
  }

  if (!is.numeric(shape2) || any(!is.finite(shape2))) {
    stop("'shape2' must be numeric and finite.")
  }

  lens <- c(length(p), length(loc), length(scale), length(shape1), length(shape2))
  L <- max(lens)

  ok_len <- function(x) length(x) %in% c(1, L)
  if (!ok_len(p) || !ok_len(loc) || !ok_len(scale) ||
      !ok_len(shape1) || !ok_len(shape2)) {
    stop("'p', 'loc', 'scale', 'shape1', and 'shape2' must have length 1 or a common length.")
  }

  p <- rep_len(p, L)
  loc <- rep_len(loc, L)
  scale <- rep_len(scale, L)
  shape1 <- rep_len(shape1, L)
  shape2 <- rep_len(shape2, L)

  xi <- shape1
  h  <- shape2

  z <- numeric(L)

  tol <- .Machine$double.eps^0.5

  # Case 1: xi != 0, h != 0
  idx11 <- abs(xi) >= tol & abs(h) >= tol
  if (any(idx11)) {
    yp <- (1 - p[idx11]^h[idx11]) / h[idx11]
    z[idx11] <- loc[idx11] +
      scale[idx11] / xi[idx11] -
      scale[idx11] / xi[idx11] * (yp^xi[idx11])
  }

  # Case 2: xi = 0, h != 0
  idx01 <- abs(xi) < tol & abs(h) >= tol
  if (any(idx01)) {
    yp <- (1 - p[idx01]^h[idx01]) / h[idx01]
    z[idx01] <- loc[idx01] - scale[idx01] * log(yp)
  }

  # Case 3: xi != 0, h = 0
  idx10 <- abs(xi) >= tol & abs(h) < tol
  if (any(idx10)) {
    yp <- -log(p[idx10])
    z[idx10] <- loc[idx10] +
      scale[idx10] / xi[idx10] -
      scale[idx10] / xi[idx10] * (yp^xi[idx10])
  }

  # Case 4: xi = 0, h = 0
  idx00 <- abs(xi) < tol & abs(h) < tol
  if (any(idx00)) {
    yp <- -log(p[idx00])
    z[idx00] <- loc[idx00] - scale[idx00] * log(yp)
  }

  z
}
