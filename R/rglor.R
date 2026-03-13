#' Random Generation from the Generalized Logistic Distribution for r-Largest Order Statistics
#'
#' Generates random samples from the generalized logistic distribution for
#' \eqn{r}-largest order statistics.
#'
#' @param n A positive integer specifying the number of observations.
#' @param r A positive integer specifying the number of order statistics
#'   for each observation.
#' @param loc A numeric value specifying the location parameter.
#' @param scale A positive numeric value specifying the scale parameter.
#' @param shape A numeric value specifying the shape parameter.
#'
#' @return A list with components:
#' \item{umat}{An \code{n x r} matrix of independent uniform random numbers.}
#' \item{wmat}{An \code{n x r} matrix of transformed uniform variables used
#' to construct decreasing order statistics.}
#' \item{rmat}{An \code{n x r} matrix of simulated \eqn{r}-largest order
#' statistics from the generalized logistic distribution.}
#'
#' @details
#' The function first generates independent uniform random variables and then
#' constructs decreasing variables through recursive transformations. These
#' are transformed using the generalized logistic quantile function
#' \code{\link{qglo}}.
#'
#'@references
#'
#' Shin, Y., & Park, J-S. (2024).
#' Generalized logistic model for r-largest order statistics with
#' hydrological application.
#' \emph{Stochastic Environmental Research and Risk Assessment}.
#' \doi{10.1007/s00477-023-02642-7}
#'
#' @export
#'
#' @examples
#' x <- rglor(10, 3, loc = 0, scale = 1, shape = 0.1)
#' x$rmat
rglor <- function(n, r, loc = 0, scale = 1, shape = 0.1) {

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

  if (!is.numeric(shape) || length(shape) != 1 || !is.finite(shape)) {
    stop("'shape' must be a finite numeric scalar.")
  }

  z <- list()

  umat <- matrix(0, nrow = n, ncol = r)
  wmat <- matrix(0, nrow = n, ncol = r)

  colnames(umat) <- paste0("u", seq_len(r))
  colnames(wmat) <- paste0("w", seq_len(r))

  i <- 1
  while (i <= n) {
    # Generate independent uniforms
    u <- stats::runif(r)

    # Construct decreasing transformed uniforms
    w <- numeric(r)
    w[1] <- u[1]

    if (r >= 2) {
      for (j in 2:r) {
        w[j] <- w[j - 1] * (u[j]^(1 / j))
      }
    }

    # Save only valid decreasing sequences
    if (all(diff(w) < 0)) {
      umat[i, ] <- u
      wmat[i, ] <- w
      i <- i + 1
    }
  }

  z$umat <- umat
  z$wmat <- wmat

  z$rmat <- matrix(
    qglo(wmat, loc = loc, scale = scale, shape = shape),
    nrow = nrow(wmat),
    ncol = ncol(wmat)
  )

  colnames(z$rmat) <- paste0("r", seq_len(r))

  invisible(z)
}
