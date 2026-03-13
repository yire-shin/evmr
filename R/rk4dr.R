#' Random Generation from the Four-Parameter Kappa Distribution for r-Largest Order Statistics
#'
#' Generates random samples from the four-parameter kappa distribution for
#' \eqn{r}-largest order statistics.
#'
#' @param n A positive integer specifying the number of observations.
#' @param r A positive integer specifying the number of order statistics
#'   for each observation.
#' @param loc A numeric value specifying the location parameter.
#' @param scale A positive numeric value specifying the scale parameter.
#' @param shape1 A numeric value specifying the first shape parameter.
#' @param shape2 A numeric value specifying the second shape parameter.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{umat}: an \code{n x r} matrix of independent uniform random numbers
#'   \item \code{wmat}: an \code{n x r} matrix of transformed uniform variables
#'   \item \code{rmat}: an \code{n x r} matrix of simulated \eqn{r}-largest order statistics
#' }
#'
#' @references
#' Shin, Y., & Park, J.-S. (2023).
#' Modeling climate extremes using the four-parameter kappa distribution
#' for r-largest order statistics.
#' \emph{Weather and Climate Extremes}.
#' \doi{10.1016/j.wace.2022.100533}
#'
#' @details
#' The function first generates independent uniform random variables and then
#' constructs decreasing transformed variables recursively using the second
#' shape parameter. These are transformed by the four-parameter kappa quantile
#' function \code{\link{qk4d}}.
#'
#' For valid generation with \eqn{r > 1}, the second shape parameter should
#' satisfy \eqn{shape2 < 1/(r-1)}.
#'
#' @export
#'
#' @examples
#' x <- rk4dr(10, 3, loc = 0, scale = 1, shape1 = 0.1, shape2 = 0.1)
#' x$rmat
rk4dr <- function(n, r, loc = 0, scale = 1, shape1 = 0.1, shape2 = 0.1) {

  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0 || n != as.integer(n)) {
    stop("'n' must be a positive integer.")
  }

  if (!is.numeric(r) || length(r) != 1 || !is.finite(r) || r <= 0 || r != as.integer(r)) {
    stop("'r' must be a positive integer.")
  }

  if (!is.numeric(loc) || length(loc) != 1 || !is.finite(loc)) {
    stop("'loc' must be a finite numeric scalar.")
  }

  if (!is.numeric(scale) || length(scale) != 1 || !is.finite(scale)) {
    stop("'scale' must be a finite numeric scalar.")
  }

  if (scale <= 0) {
    stop("'scale' must be positive.")
  }

  if (!is.numeric(shape1) || length(shape1) != 1 || !is.finite(shape1)) {
    stop("'shape1' must be a finite numeric scalar.")
  }

  if (!is.numeric(shape2) || length(shape2) != 1 || !is.finite(shape2)) {
    stop("'shape2' must be a finite numeric scalar.")
  }

  if (r > 1 && shape2 >= 1 / (r - 1)) {
    stop("'shape2' must be smaller than 1/(r-1) for valid generation.")
  }

  z <- list()

  umat <- matrix(0, nrow = n, ncol = r)
  wmat <- matrix(0, nrow = n, ncol = r)

  colnames(umat) <- paste0("u", seq_len(r))
  colnames(wmat) <- paste0("w", seq_len(r))

  i <- 1
  while (i <= n) {
    u <- stats::runif(r)

    w <- numeric(r)
    w[1] <- u[1]

    if (r >= 2) {
      for (j in 2:r) {
        expo <- 1 / (1 - (j - 1) * shape2)
        w[j] <- w[j - 1] * (u[j]^expo)
      }
    }

    if (all(diff(w) < 0)) {
      umat[i, ] <- u
      wmat[i, ] <- w
      i <- i + 1
    }
  }

  z$umat <- umat
  z$wmat <- wmat

  z$rmat <- matrix(
    qk4d(wmat, loc = loc, scale = scale, shape1 = shape1, shape2 = shape2),
    nrow = nrow(wmat),
    ncol = ncol(wmat)
  )

  colnames(z$rmat) <- paste0("r", seq_len(r))

  invisible(z)
}
