#' @name qgd
#' @aliases qgd
#' @title qunatile of gumbel
#' @param p Vector of quantiles.
#' @param loc,scale Location, scale parameters. Can be vectors, but the lengths must be appropriate.
#'
#' @return quantile function for the 1-GD distribution
#' @export
#'
#' @examples
#' \dontrun{
#' qgd(0.5, loc = 0, scale=1)
#' }
qgd<-function (p,loc = 0, scale = 1)
{
  mu  = loc
  sig = scale

  yp  = -log(p)

  z <-mu - sig * log(yp)
  z
}

