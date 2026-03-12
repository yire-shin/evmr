#' @name qglo
#' @aliases qglo
#' @title qunatile of glo
#' @param p Vector of quantiles.
#' @param loc,scale,shape Location, scale, and shape parameters. Can be vectors, but the lengths must be appropriate.
#'
#' @return quantile function for the 1-GGD distribution
#' @export
#'
#' @examples
#' \dontrun{
#' qglo(0.5, loc = 0, scale=1, shape = 0.1)
#' }
qglo<-function (p,loc = 0, scale = 1, shape = 0)
{
  mu = loc
  sig= scale
  xi = shape

  f = p
  f_inv = 1/f

  z = mu +sig/xi - sig/xi * exp(xi*log(f_inv-1))
  z
}
