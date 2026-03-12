#' @name qk4d
#' @aliases qk4d
#' @title qunatile of k4d
#' @param p Vector of quantiles.
#' @param loc,scale,shape1,shape2 Location, scale, and shape parameters. Can be vectors, but the lengths must be appropriate.
#'
#' @return quantile function for the 1-GGD distribution
#' @export
#'
#' @examples
#' \dontrun{
#' qk4d(0.5, loc = 0, scale=1, shape1 = 0.1, shape2 = 0.1)
#' }
qk4d<-function (p,loc = 0, scale = 1, shape1 = 0.1, shape2 = 0.1)
{
  mu  = loc
  sig = scale
  xi  = shape1
  h   = shape2

  yp  = (1-(p)^h) / h

  z <-mu + (sig/xi) - (sig/xi) * (yp^xi)
  z
}
