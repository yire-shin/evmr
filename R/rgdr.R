#' Random Generation from the Gumbel Distribution for r-Largest Order Statistics
#'
#' Generates random samples from the Gumbel distribution for
#' \eqn{r}-largest order statistics.
#'
#' @param n A positive integer specifying the number of observations.
#' @param r A positive integer specifying the number of order statistics
#'   for each observation.
#' @param loc A numeric value specifying the location parameter.
#' @param scale A positive numeric value specifying the scale parameter.
#'
#' @return A list with components:
#' \item{umat}{An \code{n x r} matrix of independent uniform random numbers.}
#' \item{wmat}{An \code{n x r} matrix of transformed uniform variables used
#' to construct decreasing order statistics.}
#' \item{rmat}{An \code{n x r} matrix of simulated \eqn{r}-largest order
#' statistics from the Gumbel distribution.}
#'
#' @details
#' The function first generates independent uniform random variables and then
#' constructs decreasing variables through cumulative products. These are
#' transformed using the Gumbel quantile function \code{\link{qgd}}.
#'
#' @export
#'
#' @examples
#' x <- rgdr(10, 3, loc = 0, scale = 1)
#' x$rmat
rgdr <- function(n, r, loc = 0, scale = 1) {

  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) || n <= 0 || n != as.integer(n)) {
    stop("'n' must be a positive integer.")
  }

  if (!is.numeric(r) || length(r) != 1 || !is.finite(r) || r <= 0 || r != as.integer(r)) {
    stop("'r' must be a positive integer.")
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

  z <- list()

  umat <- matrix(stats::runif(n * r), nrow = n, ncol = r)
  wmat <- umat

  if (r >= 2) {
    for (j in 2:r) {
      wmat[, j] <- wmat[, j - 1] * umat[, j]
    }
  }

  colnames(umat) <- paste0("u", seq_len(r))
  colnames(wmat) <- paste0("w", seq_len(r))

  z$umat <- umat
  z$wmat <- wmat
  z$rmat <- qgd(wmat, loc = loc, scale = scale)

  colnames(z$rmat) <- paste0("r", seq_len(r))

  invisible(z)
}
