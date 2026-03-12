#' @name rgdr
#' @aliases rgdr
#' @title random sample generation for rgd
#' @param n Number of observations
#' @param r Number of order statistics for each observation.
#' @param loc,scale Location, scale. Can be vectors, but the lengths must be appropriate.
#'
#' @return Random number generation
#' @export
#'
#' @examples
#' \dontrun{
#' rgdr(50, 5, loc = 0, scale=1)
#' }
rgdr <- function(n, r, loc = 0, scale = 1) {

  z<-list()

  umat <- matrix(0, nrow = n, ncol = r)
  wmat <- matrix(0, nrow = n, ncol = r)
  colnames(umat) <- paste0("u", 1:r)
  colnames(wmat) <- paste0("w", 1:r)

  i <- 1
  while (i <= n) {
    # uniform(0, 1) random number generation
    u <- stats::runif(r)
    # w_1 = u_1, w_2 = u_1 * u_2^(1-(2-1)*h), ..., w_r = w_(r-1) * u_r^(1-(r-1)*h)
    w <- numeric(r)
    w[1] <- u[1]
    for (j in 2:r) {
      w[j] <- w[j - 1] * (u[j])
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
  z$rmat <- qgd(wmat, loc = loc, scale = scale)

  colnames(z$rmat) <- paste0("r", 1:r)

  invisible(z)
}
