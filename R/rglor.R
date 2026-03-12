#' @name rglor
#' @aliases rglor
#' @title random sample generation for rglo
#' @param n Number of observations
#' @param r Number of order statistics for each observation.
#' @param loc,scale,shape Location, scale, and shape parameters. Can be vectors, but the lengths must be appropriate.
#'
#' @return Random number generation
#' @export
#'
#' @examples
#' \dontrun{
#' rglor(50, 5, loc = 0, scale=1, shape = 0.1)
#' }
rglor <- function(n, r, loc = 0, scale = 1, shape = 0.1) {

  z<-list()

  umat <- matrix(0, nrow = n, ncol = r)
  wmat <- matrix(0, nrow = n, ncol = r)
  colnames(umat) <- paste0("u", 1:r)
  colnames(wmat) <- paste0("w", 1:r)

  i <- 1
  while (i <= n) {
    # uniform(0, 1) random number generation
    u <- stats::runif(r)
    # w_1 = u_1, w_2 = u_1 * u_2^(r=2), ..., w_r = w_(r-1) * u_r^(r=r)
    w <- numeric(r)
    w[1] <- u[1]
    for (j in 2:r) {
      w[j] <- w[j - 1] * (u[j]^(1/j))
    }

    # save data if this condition w_1 > w_2 > ... > w_r
    if (all(diff(w) < 0)) {
      umat[i, ] <- u
      wmat[i, ] <- w
      i <- i + 1
    }
  }

  z$umat <- umat
  z$wmat <- wmat
  z$rmat <- qglo(wmat, loc = loc, scale = scale, shape = shape)

  colnames(z$rmat) <- paste0("r", 1:r)

  invisible(z)
}
